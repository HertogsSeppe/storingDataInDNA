""" 
Stack overflow, Brunão: https://stackoverflow.com/questions/75611160/
converting-any-file-to-binary-1-and-0-and-then-back-to-original-file-without-cor
"""


def file_to_binary_string(file_path):
    with open(file_path, "rb") as file:
        binary_code = file.read()
        binary_string = "".join(format(byte, "08b") for byte in binary_code)
    return binary_string


def binary_string_to_file(binary_string, file_path):
    with open(file_path, "wb") as file:
        bytes_list = [
            int(binary_string[i : i + 8], 2) for i in range(0, len(binary_string), 8)
        ]
        bytes_arr = bytearray(bytes_list)
        file.write(bytes_arr)


def DNA_srand_to_file(strand, file_path):
    with open(file_path, "w") as file:
        file.write(strand)
