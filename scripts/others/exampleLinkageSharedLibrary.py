import ctypes

def main():

    # Load the C shared library
    libfile="./libECOGEN.so"  # provide suitable path, either to make build (as done here) or cmake one ('build/libECOGEN.so')

    lib = ctypes.CDLL(libfile)

    # Call ecogen library
    # args = (c_char_p * 0)()
    lib.libecogen()

if __name__ == '__main__':
    main()
