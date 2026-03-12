import re
from collections.abc import Iterable
from typing import Optional

import pyarrow.parquet as pq
from psm_utils.io._base_classes import ReaderBase
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

RESCORING_FEATURES = [
    "Precursor.Charge",
    "RT",
    "Predicted.RT",
    "iRT",
    "Predicted.iRT",
    "Ms1.Profile.Corr",
    "Ms1.Area",
    "Ms1.Normalised",
    "IM",
    "iIM",
    "Predicted.IM",
    "Predicted.iIM",
    "Evidence",
    "Mass.Evidence",
    "Channel.Evidence",
    "Quantity.Quality",
    "Empirical.Quality",
    "Normalisation.Noise",
    "FWHM",
    "Ms1.Apex.Area",
    "Ms1.Apex.Mz.Delta",
    "Ms1.Total.Signal.Before",
    "Ms1.Total.Signal.After",
    "RT.Start",
    "RT.Stop",
]


class DIANN2ParquetReader(ReaderBase):
    def __init__(self, filename) -> None:
        """
        psm_utils Reader for DIA-NN 2.0 main report in parquet format.

        Parameters
        ----------
        file_path : str
            Path to the DIA-NN 2.0 main report.

        """
        self.filename = filename

    def __iter__(self) -> Iterable[PSM]:
        # TODO: There may have some floating-point precision issues
        # parquet stores float64
        with pq.ParquetFile(self.filename) as reader:
            for batch in reader.iter_batches():
                for row in batch.to_pylist():
                    yield self._get_psm(row)

    def _get_psm(self, psm_dict: dict) -> PSM:
        rescoring_features = {}
        for feature in RESCORING_FEATURES:
            try:
                rescoring_features[feature] = float(psm_dict[feature])
            except KeyError:
                continue

        return PSM(
            peptidoform=self._parse_peptidoform(psm_dict["Modified.Sequence"], psm_dict["Precursor.Charge"]),
            spectrum_id="NaN",  # DIA-NN 2.0 does not provide MS2 scan. Users have to align by RT.
            run=psm_dict["Run"],
            is_decoy=bool(psm_dict["Decoy"]),
            qvalue=psm_dict["Q.Value"],
            pep=float(psm_dict["PEP"]),
            score=None,  # DIA-NN 2.0 removed the "CScore" from the report
            precursor_mz=float(psm_dict["Precursor.Mz"]),
            retention_time=float(psm_dict["RT"]),
            ion_mobility=float(psm_dict["IM"]),
            protein_list=self._get_protein_list(psm_dict["Protein.Ids"]),
            source="diann2",
            rank=None,
            provenance_data=({"diann_filename": str(self.filename)}),
            rescoring_features=rescoring_features,
            metadata={},
        )

    @staticmethod
    def _parse_peptidoform(peptide: str, charge: Optional[str]) -> str:
        # Add charge
        if charge:
            peptide += f"/{int(float(charge))}"

        # Replace parentheses with square brackets and capitalize UniMod prefix
        pattern = r"\(UniMod:(\d+)\)"
        replacement = r"[UNIMOD:\1]"
        peptide = re.sub(pattern, replacement, peptide)

        # Add hyphen for N-terminal modifications
        # If [UNIMOD:n] occurs before the first amino acid, a hyphen is added before the first
        # amino acid
        if peptide[0] == "[":
            # Hyphen after the closing bracket
            peptide = peptide.replace("]", "]-", 1)

        # C-terminal modifications are currently not supported in DIA-NN

        return peptide

    @staticmethod
    def _get_protein_list(protein_ids_raw: str) -> list[str]:
        """
        DIA-NN 2.x update:
        - Protein.Ids no longer contains ';'
        - One protein group with leading "N/" (N = number of proteins)
        - Entries separated by "/"
        Example: "3/sp|Q504T8|MIDN_HUMAN/sp|P47928|ID4_HUMAN/DECOY_18802"
        """
        i = protein_ids_raw.find("/")
        if i != -1 and protein_ids_raw[0].isdigit():
            protein_ids_raw = protein_ids_raw[i + 1 :]
        return protein_ids_raw.split("/")

    @classmethod
    def from_dataframe(cls, dataframe) -> PSMList:
        return PSMList(
            ptm_list=[cls._get_peptide_spectrum_match(cls(""), entry) for entry in dataframe.to_dict(orient="records")]
        )
