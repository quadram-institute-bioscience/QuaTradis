
def create_parser(name, subparsers, function, options, help_message, description="", usage=""):
    if not description:
        description = help_message

    parser = subparsers.add_parser(name, help=help_message, description=description, usage=usage)
    options(parser)
    parser.set_defaults(func=function)

    parser.add_argument('--profile', dest='profile', action='store_true',
                        help='Turn on profiling.  Prints out cumulative time in each function to stdout and to an output file (tradis.profile).  The profile file can be read by tools such as snakeviz.')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Print out extra information about the process while it runs', default=False)

