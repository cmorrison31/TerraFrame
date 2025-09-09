from importlib import resources
import struct
import numpy as np


def main():
    # To drastically save on space, this script is used to convert an EGM2008
    # coefficient data file from text to binary. Uncertainty or error
    # coefficients are also dropped and the degree can be limited as well

    row_format = struct.Struct("<HHdd")
    file_name = 'EGM2008_to2190_TideFree'
    file_path = resources.files("TerraFrame.Data").joinpath(file_name)
    order_limit = 200
    output_file_name = f'EGM2008_to{order_limit}_TideFree.bin'

    # noinspection PyTypeChecker
    with open(file_path, 'r', encoding='utf8') as f:
        file_contents = f.readlines()

    parsed_content = []
    for line in file_contents:
        data = [float(x) for x in
                line.replace('D', 'e').split()[0:4]]

        parsed_content.append(data)

    parsed_content = np.array(parsed_content)

    with open(output_file_name, 'wb') as f:
        for row in parsed_content:
            if row[0] <= order_limit:
                f.write(row_format.pack(int(row[0]), int(row[1]),
                                        row[2], row[3]))


if __name__ == '__main__':
    main()
