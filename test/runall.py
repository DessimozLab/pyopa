import os

if __name__ == '__main__':

    darwin_input = os.path.dirname(__file__) + "/test_darwin.py"
    exec(compile(open(darwin_input, "rb").read(), darwin_input, 'exec'))
    test_input = os.path.dirname(__file__) + "/test_input.py"
    exec(compile(open(test_input, "rb").read(), test_input, 'exec'))
