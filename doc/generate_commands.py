#!/usr/bin/env python3
"""Generate the Mash command reference from the source.

Commands and options are declared in the constructors of the Command
subclasses (src/mash/Command*.cpp).  Parsing them there keeps the docs from
drifting from what 'mash <command> -h' prints.  Writes:

    doc/sphinx/commands.rst   reference page for the Sphinx manual
    doc/man/mash*.1           one man page per command, plus mash(1)

Usage:
    python3 doc/generate_commands.py [--check] [--rst-only | --man-only]
"""

import argparse
import os
import re
import struct
import sys

SRC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, 'src', 'mash')
DOC_DIR = os.path.dirname(os.path.abspath(__file__))
RST_PATH = os.path.join(DOC_DIR, 'sphinx', 'commands.rst')
MAN_DIR = os.path.join(DOC_DIR, 'man')

# Date stamped into the man page headers.  Fixed rather than today's date so
# regenerating without source changes produces no diff.
MAN_DATE = '2025-09-15'

# Category display names and their ordering, from Command::Command().
CATEGORY_DISPLAY = [
    ('', ''),
    ('Input', 'Input'),
    ('Output', 'Output'),
    ('Sketch', 'Sketching'),
    ('Window', 'Sketching (windowed)'),
    ('Reads', 'Sketching (reads)'),
    ('Alphabet', 'Sketching (alphabet)'),
]

# Argument placeholder per Option::Type, from Command::print().
TYPE_NAMES = {
    'Boolean': None,
    'Number': 'num',
    'Integer': 'int',
    'Size': 'size',
    'File': 'path',
    'String': 'text',
}

MASH_CPP = os.path.join(SRC_DIR, 'mash.cpp')

# Column at which option descriptions start in the generated option lists.
OPTION_INDENT = 13
OPTION_TEXT_WIDTH = 63

# Options that useSketchOptions() adds, in order.  The windowed ones are
# guarded by #ifdef COMMAND_FIND and so are omitted, matching a default build.
SKETCH_OPTIONS = [
    'threads', 'kmer', 'noncanonical', 'protein', 'alphabet', 'case',
    'sketchSize', 'individual', 'seed', 'warning', 'reads',
    'memory', 'minCov', 'targetCov', 'genome',
]


class Option(object):
    def __init__(self, type_, identifier, category, description,
                 default='', arg_min=0.0, arg_max=0.0):
        self.type = type_
        self.identifier = identifier
        self.category = category
        self.description = description
        self.default = default
        self.arg_min = arg_min
        self.arg_max = arg_max

    def signature(self):
        """The left-hand column, e.g. '-k <int>'."""
        text = '-' + self.identifier
        type_name = TYPE_NAMES.get(self.type)
        if type_name is not None:
            text += ' <%s>' % type_name
        return text

    def full_description(self):
        """The right-hand column, with range and default appended."""
        text = self.description
        if self.arg_min != self.arg_max:
            if self.type == 'Integer':
                text += ' (%d-%d)' % (int(self.arg_min), int(self.arg_max))
            else:
                text += ' (%s-%s)' % (format_float(self.arg_min),
                                      format_float(self.arg_max))
        if self.default != '':
            text += ' [%s]' % self.default
        return text


def format_float(value):
    """Match the way C++ ostream renders a float (e.g. 1 rather than 1.0)."""
    if value == int(value):
        return str(int(value))
    return repr(value)


class CommandDef(object):
    def __init__(self, name):
        self.name = name
        self.summary = ''
        self.description = ''
        self.argument_string = ''
        self.options = []


def strip_comments(text):
    """Remove // and /* */ comments so commented-out declarations are ignored.

    String literals are preserved, since descriptions may contain '//'.
    """
    out = []
    i = 0
    length = len(text)
    while i < length:
        char = text[i]
        if char == '"':
            j = i + 1
            while j < length:
                if text[j] == '\\':
                    j += 2
                    continue
                if text[j] == '"':
                    break
                j += 1
            out.append(text[i:j + 1])
            i = j + 1
        elif text.startswith('//', i):
            j = text.find('\n', i)
            i = length if j == -1 else j
        elif text.startswith('/*', i):
            j = text.find('*/', i + 2)
            i = length if j == -1 else j + 2
        else:
            out.append(char)
            i += 1
    return ''.join(out)


def split_arguments(text):
    """Split a C++ argument list on top-level commas."""
    args = []
    depth = 0
    current = []
    i = 0
    length = len(text)
    while i < length:
        char = text[i]
        if char == '"':
            j = i + 1
            while j < length:
                if text[j] == '\\':
                    j += 2
                    continue
                if text[j] == '"':
                    break
                j += 1
            current.append(text[i:j + 1])
            i = j + 1
            continue
        if char in '(<[':
            depth += 1
        elif char in ')>]':
            depth -= 1
        if char == ',' and depth == 0:
            args.append(''.join(current).strip())
            current = []
        else:
            current.append(char)
        i += 1
    if ''.join(current).strip():
        args.append(''.join(current).strip())
    return args


def parse_string_literal(text):
    """Evaluate a C++ string expression: adjacent literals are joined, and a
    getOption("name").identifier term resolves to that option's flag.
    """
    parts = []
    i = 0
    length = len(text)
    while i < length:
        char = text[i]
        if char == '"':
            j = i + 1
            chunk = []
            while j < length:
                if text[j] == '\\':
                    chunk.append(text[j + 1])
                    j += 2
                    continue
                if text[j] == '"':
                    break
                chunk.append(text[j])
                j += 1
            parts.append(('literal', ''.join(chunk)))
            i = j + 1
            continue
        match = re.match(r'getOption\(\s*"([^"]+)"\s*\)\s*\.\s*identifier', text[i:])
        if match:
            parts.append(('option', match.group(1)))
            i += match.end()
            continue
        i += 1
    return parts


def parse_number(text):
    """Parse a numeric literal, narrowed to float as the C++ stores it.

    Option bounds are floats, so 0xFFFFFFFF rounds to 4294967296 -- which is
    what the binary prints.
    """
    text = text.strip().rstrip('.')
    if not text:
        return 0.0
    try:
        value = float(int(text, 0))
    except ValueError:
        try:
            value = float(text)
        except ValueError:
            return 0.0
    return struct.unpack('f', struct.pack('f', value))[0]


def parse_registered_classes(path):
    """Return the Command classes mash.cpp registers in a default build.

    Commands guarded by an #ifdef (currently CommandContain and CommandFind)
    are not compiled in and so are left out of the reference.
    """
    with open(path) as handle:
        source = strip_comments(handle.read())

    registered = []
    depth = 0
    for line in source.split('\n'):
        stripped = line.strip()
        if re.match(r'#\s*if', stripped):
            depth += 1
            continue
        if re.match(r'#\s*endif', stripped):
            depth = max(0, depth - 1)
            continue
        if re.match(r'#\s*(else|elif)', stripped):
            continue
        if depth > 0:
            continue
        match = re.search(r'addCommand\s*\(\s*new\s+(?:mash::)?(\w+)', stripped)
        if match:
            registered.append(match.group(1))
    return registered


def parse_available_options(command_cpp):
    """Read the shared options declared in Command::Command()."""
    available = {}
    for match in re.finditer(
            r'addAvailableOption\s*\(\s*"([^"]+)"\s*,\s*Option\s*\((.*?)\)\s*\)\s*;',
            command_cpp, re.DOTALL):
        name = match.group(1)
        option = build_option(match.group(2), {})
        if option is not None:
            available[name] = option
    return available


def build_option(args_text, options_by_name):
    """Construct an Option from the arguments of an Option(...) constructor."""
    args = split_arguments(args_text)
    if len(args) < 4:
        return None

    type_ = args[0].split('::')[-1].strip()

    identifier_parts = parse_string_literal(args[1])
    if not identifier_parts:
        return None
    identifier = identifier_parts[0][1]

    category_parts = parse_string_literal(args[2])
    category = category_parts[0][1] if category_parts else ''

    description = ''
    for kind, value in parse_string_literal(args[3]):
        if kind == 'literal':
            description += value
        elif value in options_by_name:
            description += options_by_name[value].identifier

    default = ''
    if len(args) > 4:
        default_parts = parse_string_literal(args[4])
        if default_parts:
            default = default_parts[0][1]

    arg_min = parse_number(args[5]) if len(args) > 5 else 0.0
    arg_max = parse_number(args[6]) if len(args) > 6 else 0.0

    return Option(type_, identifier, category, description, default, arg_min, arg_max)


def parse_command(path, available):
    """Parse one Command subclass constructor into a CommandDef."""
    with open(path) as handle:
        source = strip_comments(handle.read())

    name_match = re.search(r'\bname\s*=\s*"([^"]+)"\s*;', source)
    if not name_match:
        return None

    command = CommandDef(name_match.group(1))

    summary_match = re.search(r'\bsummary\s*=\s*(.*?);', source, re.DOTALL)
    if summary_match:
        command.summary = ''.join(
            value for kind, value in parse_string_literal(summary_match.group(1)))

    description_match = re.search(r'\bdescription\s*=\s*(.*?);', source, re.DOTALL)
    if description_match:
        command.description = ''.join(
            value for kind, value in parse_string_literal(description_match.group(1)))

    argument_match = re.search(r'\bargumentString\s*=\s*(.*?);', source, re.DOTALL)
    if argument_match:
        command.argument_string = ''.join(
            value for kind, value in parse_string_literal(argument_match.group(1)))

    # Walk the constructor in order, so options appear as the binary lists them.
    constructor = source[name_match.start():]
    end = constructor.find('\n}')
    if end != -1:
        constructor = constructor[:end]

    options_by_name = {}
    pattern = re.compile(
        r'(useSketchOptions\s*\(\s*\)|useOption\s*\(\s*"([^"]+)"\s*\)|'
        r'addOption\s*\(\s*"([^"]+)"\s*,\s*Option\s*\((.*?)\)\s*\)\s*;)',
        re.DOTALL)

    for match in pattern.finditer(constructor):
        if match.group(1).startswith('useSketchOptions'):
            for option_name in SKETCH_OPTIONS:
                if option_name in available:
                    add_option(command, options_by_name, option_name,
                               available[option_name])
        elif match.group(2):
            option_name = match.group(2)
            if option_name in available:
                add_option(command, options_by_name, option_name,
                           available[option_name])
        else:
            option = build_option(match.group(4), options_by_name)
            if option is not None:
                add_option(command, options_by_name, match.group(3), option)

    return command


def add_option(command, options_by_name, name, option):
    if name in options_by_name:
        return
    options_by_name[name] = option
    command.options.append(option)


def wrap(text, width, indent=''):
    """Wrap text to a width, without breaking words."""
    words = text.split()
    lines = []
    current = ''
    for word in words:
        candidate = word if not current else current + ' ' + word
        if len(candidate) > width and current:
            lines.append(indent + current)
            current = word
        else:
            current = candidate
    if current:
        lines.append(indent + current)
    return lines


def render(commands):
    out = []
    out.append('.. This file is generated by doc/generate_commands_rst.py.')
    out.append('   Do not edit it directly; edit the option definitions in')
    out.append('   src/mash/Command*.cpp and regenerate.')
    out.append('')
    out.append('Command reference')
    out.append('=================')
    out.append('')
    out.append('Every :code:`mash` operation is a subcommand of the single')
    out.append(':code:`mash` executable::')
    out.append('')
    out.append('  mash <command> [options] [arguments ...]')
    out.append('')
    out.append('The same reference is available at the terminal: :code:`mash`')
    out.append('lists the commands, and :code:`mash <command> -h` describes one')
    out.append('command and its options.')
    out.append('')

    # Summary table of all commands.
    name_width = max(len(c.name) for c in commands)
    name_width = max(name_width, len('Command'))
    summary_width = max(len(c.summary) for c in commands)
    summary_width = max(summary_width, len('Description'))
    rule = '=' * name_width + ' ' + '=' * summary_width

    out.append(rule)
    out.append('Command'.ljust(name_width) + ' ' + 'Description')
    out.append(rule)
    for command in commands:
        out.append(command.name.ljust(name_width) + ' ' + command.summary)
    out.append(rule)
    out.append('')

    for command in commands:
        title = 'mash ' + command.name
        out.append('')
        out.append(title)
        out.append('-' * len(title))
        out.append('')
        out.append('Usage::')
        out.append('')
        usage = '  mash %s [options]' % command.name
        if command.argument_string:
            usage += ' ' + command.argument_string
        out.append(usage)
        out.append('')
        out.extend(wrap(command.description, 76))
        out.append('')

        if not command.options:
            continue

        out.append('Options:')
        out.append('')

        current_category = None
        for category, display in CATEGORY_DISPLAY:
            options = [o for o in command.options if o.category == category]
            if not options:
                continue
            if display:
                out.append('')
                out.append('*%s*' % display)
                out.append('')
            for option in options:
                # An RST option list: the description is indented to a
                # fixed column so Sphinx renders a two-column table.
                signature = option.signature()
                lines = wrap(option.full_description(),
                             OPTION_TEXT_WIDTH, ' ' * OPTION_INDENT)
                first = lines[0].strip() if lines else ''
                if len(signature) + 2 > OPTION_INDENT:
                    out.append(signature)
                    out.append(' ' * OPTION_INDENT + first)
                else:
                    out.append(signature.ljust(OPTION_INDENT) + first)
                out.extend(lines[1:])
                out.append('')

    return '\n'.join(out).rstrip() + '\n'


def man_escape(text):
    r"""Escape text for nroff.

    Backslashes become \e, and hyphens are written \- so they render as
    literal minus signs rather than being treated as hyphenation points.
    """
    return text.replace('\\', r'\e').replace('-', r'\-')


def man_header(title, date):
    """The page preamble, matching what Asciidoctor produced previously."""
    return [
        '\'\\" t',
        r'.\"     Title: %s' % title,
        r'.\"    Author: [see the "AUTHOR(S)" section]',
        r'.\" Generator: doc/generate_commands.py',
        r'.\"      Date: %s' % date,
        '.\\"    Manual: \\ \\&',
        '.\\"    Source: \\ \\&',
        r'.\"  Language: English',
        r'.\"',
        '.TH "%s" "1" "%s" "\\ \\&" "\\ \\&"' % (man_escape(title.upper()), date),
        r'.ie \n(.g .ds Aq \(aq',
        ".el       .ds Aq '",
        r'.ss \n[.ss] 0',
        '.nh',
        '.ad l',
        '.de URL',
        r'\fI\\$2\fP <\\$1>\\$3',
        '..',
        '.als MTO URL',
        '.if \\n[.g] \\{\\',
        '.  mso www.tmac',
        '.  am URL',
        '.    ad l',
        '.  .',
        '.  am MTO',
        '.    ad l',
        '.  .',
        '.  LINKSTYLE blue R < >',
        '.\\}',
    ]


def render_man_command(command, date):
    """Render the man page for a single 'mash <command>'."""
    page = 'mash-' + command.name
    out = man_header(page, date)

    summary = command.summary.rstrip('.')
    if summary:
        summary = summary[0].lower() + summary[1:]

    out.append('.SH "NAME"')
    out.append('%s \\- %s' % (man_escape(page), man_escape(summary)))

    out.append('.SH "SYNOPSIS"')
    out.append('.sp')
    synopsis = r'\fBmash %s\fP [options]' % man_escape(command.name)
    if command.argument_string:
        synopsis += ' ' + man_escape(command.argument_string)
    out.append(synopsis)

    out.append('.SH "DESCRIPTION"')
    out.append('.sp')
    out.extend(wrap(man_escape(command.description), 72))

    if command.options:
        out.append('.SH "OPTIONS"')
        for category, display in CATEGORY_DISPLAY:
            options = [o for o in command.options if o.category == category]
            if not options:
                continue
            if display:
                out.append('.SS "%s"' % display)
            for option in options:
                out.append('.sp')
                signature = r'\fB\-%s\fP' % man_escape(option.identifier)
                type_name = TYPE_NAMES.get(option.type)
                if type_name is not None:
                    signature += ' <%s>' % type_name
                out.append(signature)
                out.append('.RS 4')
                out.extend(wrap(man_escape(option.full_description()), 72))
                out.append('.RE')

    out.append('.SH "SEE ALSO"')
    out.append('.sp')
    out.append('mash(1)')
    return '\n'.join(out) + '\n'


def render_man_main(commands, date):
    """Render the top-level mash(1) page, listing every command."""
    out = man_header('mash', date)

    out.append('.SH "NAME"')
    out.append('mash \\- fast genome and metagenome distance estimation '
               'using MinHash')

    out.append('.SH "SYNOPSIS"')
    out.append('.sp')
    out.append(r'\fBmash\fP <command> [options] [arguments ...]')

    out.append('.SH "DESCRIPTION"')
    out.append('.sp')
    out.extend(wrap(man_escape(
        'mash is the main executable for the Mash software. The actual '
        "functionality is provided by the subtools ('commands'):"), 72))

    out.append('.SS "Commands"')
    for command in commands:
        out.append('.sp')
        out.append(r'\fB%s\fP' % man_escape(command.name))
        out.append('.RS 4')
        out.extend(wrap(man_escape(command.summary), 72))
        out.append('.RE')

    out.append('.SH "SEE ALSO"')
    out.append('.sp')
    out.append(',\n'.join('mash%s(1)' % man_escape('-' + c.name)
                          for c in commands))
    return '\n'.join(out) + '\n'


def parse_commands():
    """Parse every registered Command subclass into a CommandDef."""
    with open(os.path.join(SRC_DIR, 'Command.cpp')) as handle:
        command_cpp = strip_comments(handle.read())
    available = parse_available_options(command_cpp)
    registered = parse_registered_classes(MASH_CPP)

    commands = []
    for entry in sorted(os.listdir(SRC_DIR)):
        if not re.match(r'Command[A-Z]\w*\.cpp$', entry):
            continue
        if entry == 'CommandList.cpp':
            continue
        if entry[:-len('.cpp')] not in registered:
            continue
        command = parse_command(os.path.join(SRC_DIR, entry), available)
        if command is not None:
            commands.append(command)

    commands.sort(key=lambda c: c.name)
    return commands


def write_or_check(path, text, check, stale):
    """Write text to path, or record it in 'stale' if --check and it differs."""
    if check:
        existing = None
        if os.path.exists(path):
            with open(path) as handle:
                existing = handle.read()
        if existing != text:
            stale.append(path)
        return False
    directory = os.path.dirname(path)
    if directory and not os.path.isdir(directory):
        os.makedirs(directory)
    with open(path, 'w') as handle:
        handle.write(text)
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Generate the Mash command reference from the source.')
    parser.add_argument('--check', action='store_true',
                        help='exit non-zero if any generated file is out of date')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--rst-only', action='store_true',
                       help='only generate the Sphinx command reference')
    group.add_argument('--man-only', action='store_true',
                       help='only generate the man pages')
    args = parser.parse_args()

    commands = parse_commands()
    if not commands:
        sys.stderr.write('error: no commands found in %s\n' % SRC_DIR)
        return 1

    outputs = []
    if not args.man_only:
        outputs.append((RST_PATH, render(commands)))
    if not args.rst_only:
        outputs.append((os.path.join(MAN_DIR, 'mash.1'),
                        render_man_main(commands, MAN_DATE)))
        for command in commands:
            outputs.append((os.path.join(MAN_DIR, 'mash-%s.1' % command.name),
                            render_man_command(command, MAN_DATE)))

    stale = []
    written = 0
    for path, text in outputs:
        if write_or_check(path, text, args.check, stale):
            written += 1

    if args.check:
        if stale:
            sys.stderr.write('error: out of date, run %s:\n'
                             % os.path.basename(__file__))
            for path in stale:
                sys.stderr.write('  %s\n' % os.path.relpath(path, DOC_DIR))
            return 1
        sys.stdout.write('up to date: %d files, %d commands\n'
                         % (len(outputs), len(commands)))
        return 0

    sys.stdout.write('wrote %d files for %d commands\n' % (written, len(commands)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
