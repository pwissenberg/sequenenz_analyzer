"""Class for the calcualtion of the spots."""
import pandas as pd


class SequenzAnalyzer:
    """Represents the Data and offers differen functions to calculate the spots."""

    def __init__(
        self,
        amino_acids: str,
        dssp_3_pred: str,
        dssp_8_pred: str,
        conserv_pred: str,
        membran_tmbed_pred: str,
        # seth_disorder_pred: str,
        bind_metall_pred: str,
        bind_nucleic_pred: str,
        bind_small_pred: str,
        protein_name: str,
    ):
        """
        Initalize the sequenz_analyzer class

        Args:
            amino_acdis(str): a string representing the amino acid string.
            dssp_3_pred(str): a string representing the DSSP3-Prediction for each amino acid.
            coonserv_pred(str): a string representing the conservation score for each amino acid.
            membran_tmbed_pred(str): a string representing the membrane tm-bed for each amino acid.
            seth_disorder_pred(str): a string representing the seth disorder for each amino acid.
            bind_metall_pred(str): a string representing the binding metal for each amino acid.
            bind_nucleic_pred(str): a string representing the binding nucleic for each amino acid.
            bind_small_nucleic_pred(str): a string representing the binding.
            protein_name(str): a string representing the name of the protein.

        Raises:
            ValueError: if one of the propeties do not have the same lenght!

        Example:
            >>> SequenzAnalyzer(....)
            ...
        """
        if (
            len(amino_acids) != len(dssp_3_pred)
            or len(amino_acids) != len(conserv_pred)
            or len(dssp_3_pred) != len(dssp_8_pred)
            or len(membran_tmbed_pred) != len(amino_acids)
            # or len(seth_disorder_pred) != len(amino_acids)
            or len(bind_metall_pred) != len(amino_acids)
            or len(bind_metall_pred) != len(bind_nucleic_pred)
            or len(bind_nucleic_pred) != len(bind_small_pred)
        ):
            wrong_seq = self.find_sequence_length_order(
                [
                    amino_acids,
                    dssp_3_pred,
                    dssp_8_pred,
                    conserv_pred,
                    membran_tmbed_pred,
                    # seth_disorder_pred,
                    bind_metall_pred,
                    bind_nucleic_pred,
                    bind_small_pred,
                ]
            )
            raise ValueError(
                f"The sequenzes do not have the same length! This seq. creates the error: {wrong_seq}"
            )
        self.protein_name = protein_name
        self.data = [
            [char for char in amino_acids],
            [char for char in dssp_3_pred],
            [char for char in dssp_8_pred],
            [char for char in conserv_pred],
            [char for char in membran_tmbed_pred],
            # [char for char in seth_disorder_pred],
            [char for char in bind_metall_pred],
            [char for char in bind_nucleic_pred],
            [char for char in bind_small_pred],
        ]
        self.protein_name: str
        self.occurence_dict = {}
        self.positions_dict = {}

        self.aas = self.data[0]
        self.dssp_3 = self.data[1]
        self.dssp_8 = self.data[2]
        self.conserv = self.data[3]
        self.membran_tmbed = self.data[4]
        # self.seth_disorder = self.data[5]
        self.bind_metall = self.data[5]
        self.bind_nucleic = self.data[6]
        self.bind_small = self.data[7]

    def calculate_spots(self):
        """
        Goes through all sequenzes and adds the comination to the dictionary.
        If the combiantion already exists it increases the counter.

        Raises:
            No Exceptions!

        Example:
            >>> calculate_spots(...)
            ...
        """
        for index, char in enumerate(self.aas):
            combination = (
                f"{self.aas[index]}{self.dssp_3[index]}{self.dssp_8[index]}"
                + f"{self.conserv[index]}{self.membran_tmbed[index]}"
                + f"{self.bind_metall[index]}{self.bind_nucleic[index]}"
                + f"{self.bind_small[index]}"
            )
            if combination in self.occurence_dict:
                self.occurence_dict[combination] += 1
                self.positions_dict[combination].append(index)
            else:
                self.occurence_dict[combination] = 1
                self.positions_dict[combination] = [index]

    def write_spots_file(self):
        """
        Saves the spots_ with count to a csv file.

        Raises:
            No Exceptions!
        """
        file_name = self.protein_name.replace(" ", "_")
        file_name += "_occ.csv"
        with open("occurences/" + file_name, "w", encoding="utf-8") as f:
            for key, value in self.occurence_dict.items():
                f.write(key + "," + str(value) + "\n")
        file_name = self.protein_name.replace(" ", "_")
        file_name += "_pos.csv"
        with open("positions/" + file_name, "w", encoding="utf-8") as f:
            for key, value in self.positions_dict.items():
                f.write(key + "," + str(value) + "\n")

    def find_sequence_length_order(self, sequences: list[str]) -> str:
        """
        Returns the sequence length order of the protein.

        Args:
            sequences(list(str)): which has to be checked.

        Returns:
            Which string has not the same length.

        Raises:
            No Exceptions!

        Example:
            >>> sequenz_analyzer.find_sequence_lenght_order([..])
            ..

        """
        for index, seq in enumerate(sequences):
            if len(sequences[0]) != len(seq):
                return f"{index}: {seq}"
        return ""


def main():
    # Read in the dataframe
    df = pd.read_csv("merged_df.tsv", sep="\t")

    for index, row in df.iterrows():
        # Create the object
        seq_analzyer = SequenzAnalyzer(
            row["seqs"],
            row["dssp3_pred"],
            row["dssp8_pred"],
            row["conservation_pred"],
            row["membrane_tmbed"],
            # row["seth_disorder_pred"],
            row["binding_bindEmbed_metal_pred"],
            row["binding_bindEmbed_nucleic_pred"],
            row["binding_bindEmbed_small_pred"],
            protein_name=row["loc"] + "_" + row["ids"].split("|")[1],
        )
        seq_analzyer.calculate_spots()
        sorted_items = sorted(
            seq_analzyer.occurence_dict.items(), key=lambda x: x[1], reverse=True
        )
        for key, value in sorted_items:
            print(f"{key}: {value}")
        seq_analzyer.write_spots_file()
        print(row["ids"].split("|")[1])


if __name__ == "__main__":
    main()
