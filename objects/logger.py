import datetime
import os


class Logger:

    def __init__(self, name: str, silent: bool, sub_prog: str = None, args=None):
        """
        A constructor.
        :param name: The analysis name
        :param silent: The silent flag
        :param sub_prog: The sub-program name
        :param args: The input arguments
        """
        self.name = name
        self.silent = silent
        self.err_count = 0
        self.warning_count = 0
        self.logfile = f'{name}.log'
        self.args = args
        if not sub_prog:
            sub_prog = ''
        self.subprog = sub_prog
        self.start()
        self.logfile = os.path.abspath(self.logfile)

    def start(self):
        """
        Adds starting information.
        :return: None.
        """
        line = f'Starting {self.subprog}'
        self.info(line)
        if self.args:
            input_args = 'Input arguments: '
            args = vars(self.args)
            for var in args:
                if args[var]:
                    arg_val = args[var]
                    if isinstance(arg_val, list):
                        arg_val = f'{arg_val[0]} and {len(arg_val) - 1} more files'
                    input_args = f'{input_args}{var}: {arg_val}\t'
            line = f'{self.subprog}\t{self.name}\tINFO\t{datetime.datetime.now()}\t{input_args}'
            self.write_to_logfile(line)

    def info(self, msg: str, command=False):
        """
        Adds information to the log file.
        :param msg: The message to write.
        :param command: A command flag. If True - do not print to screen, to prevent an excessive output.
        :return:
        """
        if (not self.silent) and (not command):
            print(msg)

        line = f'{datetime.datetime.now()}\t{msg}'
        line = f'{self.subprog}\t{self.name}\tINFO\t{line}'
        self.write_to_logfile(line)

    def error(self, msg: str):
        """
        Writes an error to the log.
        :param msg: The error message to write.
        :return: None
        """
        self.err_count += 1

        print('Error:', msg)  # Write even if silence

        msg = msg.replace('\n', ' ')  # Replacing the newline character
        line = f'ERROR\t{datetime.datetime.now()}\t{msg}'

        self.write_to_logfile(f'{self.subprog}\t{self.name}\t{line}')

    def warning(self, msg: str):
        """
        Writes a warning to the log.
        :param msg: The warning message to write.
        :return: None
        """
        self.warning_count += 1

        if not self.silent:
            print('Warning:', msg)

        msg = msg.replace('\n', ' ')  # Replacing the newline character
        line = f'WARNING\t{datetime.datetime.now()}\t{msg}'

        self.write_to_logfile(f'{self.subprog}\t{self.name}\t{line}')

    def done(self):
        """
        Writes the closing comments.
        :return: None
        """
        if (self.err_count == 0) and (self.warning_count == 0):
            error_str = 'Finished successfully'
        else:
            error_str = f'{self.err_count} errors and {self.warning_count} warnings occurred'
        line = f'Done. {error_str}. Log file: {os.path.relpath(self.logfile)}\n'
        self.info(line)

    def write_to_logfile(self, line: str):
        """
         Writes a message to the logfile.
        :param line: The line to write to the logfile.
        :return:
        """
        line = f'{line}\n'  # Adding a 'new line' character
        try:
            with open(self.logfile, 'a') as f:
                f.write(line)
        except IOError as e:
            print(f'Error: could not write to the log file {self.logfile}. Error code: {e}')

    def move_log(self, directory: str):
        """
        Moves the log file to a different directory.
        :param directory: The new directory name.
        :return: None.
        """
        new_path = f'{directory}/{os.path.basename(self.logfile)}'
        old_path = self.logfile
        if os.path.abspath(new_path) != os.path.abspath(old_path):  # The file name are different
            try:
                if os.path.isfile(new_path):  # The file exists
                    with open(new_path, 'a') as nf:  # Append to the new file
                        with open(old_path) as of:
                            for line in of:
                                nf.write(line)
                    os.remove(old_path)

                else:
                    os.rename(self.logfile, new_path)  # Replace the file name

                self.logfile = os.path.abspath(new_path)
            except IOError as e:
                print(f'Error: could not move the log file {self.logfile} to the directory {dir}: {e}')

    def update_errors(self, count: int):
        """
        Updates the errors count
        :param count: The count to add to the errors count
        :return: None
        """
        if count > 0:
            self.err_count += count

    def update_warnings(self, count: int):
        """
        Updates the warnings count
        :param count: The count to add to the warnings count
        :return: None
        """
        if count > 0:
            self.warning_count += count
