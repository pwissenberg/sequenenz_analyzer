"""
This class helps to find the best Amino acid' properties combination to separate in endo and exo.
"""
import os
from enum import Enum
from itertools import combinations

import pandas as pd

pd.options.display.max_rows = 100


class Embedding(Enum):
    """Represents the different positions of the embeddings."""

    AMINO_ACID = 0
    DSSP_3 = 1
    DSSP_8 = 2
    CONSVERATION_SCORE = 3
    MEMBRANE_TMBED = 4
    BIND_METALL = 5
    BIND_NUCLEIC = 6
    BIND_SM_MOL = 7
    POSITION = 8


class AminoAcidFinder:
    """Contain all the different properties."""

    def __init__(self) -> None:
        """
        Constructs an AminoAcidFinder Object

        Args:


        Returns:


        Raises:
            No Exceptions!

        Example:
            >>>


        """

    def calc_possible_permuations(self) -> None:
        """
        Give back all the different permutations of the combinations.
        Args:

        Returns:


        Raises:
            No Exceptions!

        Example:
            >>>


        """
        pass

    def calc_occurence(
        self, permutation: str, df: pd.DataFrame, prot_id: str, loc: str
    ) -> pd.DataFrame:
        """
        Calculate the occurence with the used permutation.

        Args:
            permutation(str): defines index should be used.


        Returns:


        Raises:
            No Exceptions!

        Example:
            >>>

        """
        word_dict = df.to_dict(orient="records")
        word_set = {}
        result_df = pd.DataFrame(
            columns=["perm", "prot_id", "amino_word", "count", "loc"]
        )

        for entry in word_dict:
            amino_acid_word = entry["word"]
            perm_key = self.get_index_amino_acid(amino_acid_word, permutation)
            if perm_key in word_set:
                word_set[perm_key].append(amino_acid_word)
            else:
                word_set[perm_key] = [amino_acid_word]

        for elem in word_set.keys():
            occurence = 0
            for amino in word_set[elem]:
                n = df.loc[df["word"] == amino, "n"]
                occurence += int(n)
            new_row = pd.DataFrame(
                {
                    "perm": permutation,
                    "prot_id": prot_id,
                    "amino_word": elem,
                    "count": occurence,
                    "loc": loc,
                },
                index=[0],
            )
            result_df = pd.concat([result_df, new_row], ignore_index=True)

        return result_df

    def get_index_amino_acid(self, amino_acid: str, permutation: str):
        return "".join([amino_acid[int(i)] for i in permutation])

    def evaluate_best_combinations(self, df: pd.DataFrame):
        """
        Evaluate the best combination of properties to separate in endo and exo.
        Thereby, it calculates the occurence of the combinaton and compares then in both groups.

        Args:


        Returns:


        Raises:
            No Exceptions!

        Example:
            >>>


        """
        # Order the rows
        pass

    def change_embedding_values(self):
        """
        Change embedding with new value heuristik for the value.

        Args:


        Returns:


        Raises:
            No Exceptions!

        Example:
            >>>


        """
        pass

    def translate_position(self):
        pass

    def read_in_permuation_files(self):
        """"""
        pass


def generate_substrings(embedding_str: str, k: int) -> set:
    """
    Calculate all the possible combinations of the embeddings to describe an amino acid.

    Args:
        embedding_str(str): Our embeddings represented as one string.
        removed_embed(str): String we want do leave out of the observation.

    Returns:
        list(str: Contains all the different

    Raises:
        No Exceptions!

    Example:
        >>>
    """
    length = len(embedding_str)
    return set(["".join(i) for i in combinations(embedding_str, length - k)])


def calc_all_permuations():
    all_permutatiomns = set()
    for i in range(0, 8):
        runner_set = generate_substrings("01234567", i)
        all_permutatiomns = all_permutatiomns | runner_set
    return all_permutatiomns


def main():
    # Read in all files
    os.chdir("/Users/paulwissenberg/Coding/sequenenz_analyzer/occurences")
    files = os.listdir()

    # all permutations
    perm_set = calc_all_permuations()
    finder = AminoAcidFinder()
    for perm in perm_set:
        # Create df
        permuation_df = pd.DataFrame(
            columns=["perm", "prot_id", "amino_word", "count", "loc"]
        )
        for index, f in enumerate(files):
            # Read in the file
            prot_df = pd.read_csv(f, sep=",", header=None, names=["word", "n"])

            prot_id = f.split("_")[1]
            loc = f.split("_")[0]

            prot_df = finder.calc_occurence(perm, prot_df, prot_id, loc)

            # Append at the df
            permuation_df = pd.concat([permuation_df, prot_df], ignore_index=True)

        # Save df
        file_path = (
            f"/Users/paulwissenberg/Coding/sequenenz_analyzer/permuation_dfs/{perm}.csv"
        )
        permuation_df.to_csv(file_path, sep=",", index=False)

        # Evaluate
        print(permuation_df)
        amino_counts = (
            permuation_df.groupby("amino_word")["prot_id"].nunique().reset_index()
        )
        amino_counts = amino_counts.sort_values(by="prot_id", ascending=False)
        endo_amino_counts = amino_counts[amino_counts["loc"] == "endo"]
        exo_amino_counts = amino_counts[amino_counts["loc"] == "exo"]

        print(amino_counts)

        print(endo_amino_counts)
        print(exo_amino_counts)

        print("______")
        print(endo_amino_counts.sort_values(by="prot_id", ascending=False))
        print(exo_amino_counts.sort_values(by="prot_id", ascending=False))

        # Make an evaluation entry

        break


if __name__ == "__main__":
    main()
