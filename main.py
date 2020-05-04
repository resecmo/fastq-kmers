import argparse
import csv
from os.path import isfile
from sys import stderr
import re
from typing import List, Tuple, NoReturn, Iterable


def csv_writer(data: Iterable, path: str) -> NoReturn:
    with open(path, "w", newline='') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(["K-mer", "Frequency"])
        for line in data:
            writer.writerow(line)


def sequences_from_fastq(filename: str) -> List[str]:
    # regexp = re.compile(r"@([A-Za-z0-9_.: -/]+)\n([A-Za-z\n.~]+)\n\+\1{0,1}\n[!-~\n]+")
    regexp = re.compile(r"@[A-Za-z0-9_.: /\-]+\n([A-Za-z\n.~]+)\n\+")
    result = []
    with open(filename, 'r') as file:
        for seq in re.findall(regexp, file.read()):
            result.append("".join(seq.split()))
    return result


def spectrum_from_sequences(seqs: Iterable[str], k: int) -> List[Tuple[str, int]]:
    spectrum = {}
    for sequence in seqs:
        for i in range(len(sequence) - k + 1):
            if sequence[i:i+k] in spectrum:
                spectrum[sequence[i:i+k]] += 1
            else:
                spectrum[sequence[i:i+k]] = 1
    return sorted(spectrum.items(), key=lambda x: x[1], reverse=True)  # sorted by frequency


def spectrum_from_fastq(input_file: str, output_file: str, k: int, force: bool = False) -> NoReturn:
    if not force and isfile(output_file):
        print("Output file already exists\n", file=stderr)
        return

    sequences = sequences_from_fastq(input_file)
    spectrum = spectrum_from_sequences(sequences, k)
    csv_writer(spectrum, output_file)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="FASTQ file from which sequences are read")
    parser.add_argument("output_file", help="Output csv file")
    parser.add_argument("k", type=int, help="Length of k-mers")
    parser.add_argument("-f", "--force", action="store_true", help="Allows file overwriting")

    args = parser.parse_args()

    spectrum_from_fastq(args.input_file, args.output_file, args.k, args.force)
