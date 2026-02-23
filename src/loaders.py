from __future__ import annotations

import functools
from typing import List, NamedTuple, Optional
from pathlib import Path

import os
import re

import anndata as ad
import pandas as pd
import scanpy as sc
import numpy as np

from attrs import frozen

from abc import ABC, abstractmethod

# Opt-in to the future behavior
pd.set_option("future.no_silent_downcasting", True)


class BaseDataCleaner(ABC):
    """Clean Tabular Data. Class to remove typos, preprocess data columns.

    Attributes
    ----------
    dataset_name : str

    Methods
    -------
    transform: Transform the input data.
    """

    @abstractmethod
    def transform(self, x):
        """Transform the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos
        """
        pass


class MOSAICLightClinicalDataCleaner(BaseDataCleaner):
    """Add sample_id as index and remove columns with only nan values.

    Useful only for the original clinical data from MOSAIC.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos
        """
        # Add sample_id as index
        x.set_index("sample_id", inplace=True)

        # Remove columns with only nan values
        x.dropna(axis=1, how="all", inplace=True)

        return x


class MOSAICClinicalDataCleaner(BaseDataCleaner):
    """Clean processed clinical data from MOSAIC.

    Remove typos, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos
        """
        # Add sample_id as index
        x.set_index("sample_id", inplace=True)

        # Replace Unknown values by nan
        x.replace("Unknown", np.nan, inplace=True)

        # Remove columns with only nan values
        x.dropna(axis=1, how="all", inplace=True)

        return x


class MOSAICKeyEventsDataCleaner(BaseDataCleaner):
    """Clean Key Events data file from MOSAIC.

    Remove typos, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos
        """
        # Add Patient ID as index
        x.set_index("patient_id", inplace=True)

        # Replace Unknown values by nan
        x.replace("Unknown", np.nan, inplace=True)

        # Remove columns with only nan values
        x.dropna(axis=1, how="all", inplace=True)

        return x


class MOSAICBulkRNADataCleaner(BaseDataCleaner):
    """Clean Key Events data file from MOSAIC.

    Remove typos, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos
        """
        x.rename(columns={"gene": "EnsemblID"}, inplace=True)
        x.set_index("EnsemblID", inplace=True)

        # If weird ensembl IDs exist, drop them
        # Appears in the deseq normalized counts for instance
        x = x[~x.index.str.contains("g@-ext")]  # corresponds to
        # immunoglobulin super-loci, mainly used to find fusion genes
        # and therefore not necessary for DGEA

        return x


class MOSAICWESDataCleaner(BaseDataCleaner):
    """Clean WES data file from MOSAIC.

    Remove typos, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos
        """
        x = x.dropna(subset=["gene_name"])
        x = x.set_index("gene_name")
        x = x.drop(columns=["ensembl_gene_id"])
        x = x.T  # Transpose the data to have patient id in index
        x.index.name = "sample_id"

        return x


class MOSAICHEDataCleaner(BaseDataCleaner):
    """Clean HE data file from MOSAIC.

    Add patient ID, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos, parsed dates
        """
        x.loc[:, "Subject Id"] = x["path"].apply(lambda v: v.name[:9])
        x.loc[:, "patient id"] = x["Subject Id"].apply(lambda v: v[:-1])
        x = x.set_index("Subject Id")
        return x


class MOSAICHEFeaturesDataCleaner(BaseDataCleaner):
    """Clean HE features data file from MOSAIC.

    Add patient ID, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos, parsed dates
        """
        x.loc[:, "Subject Id"] = x["path"].apply(lambda v: v.name.replace("features_", "")[:9])
        x.loc[:, "patient id"] = x["Subject Id"].apply(lambda v: v[:-1])
        x = x.set_index("Subject Id")
        return x


class BruceLightMetadataDataCleaner(BaseDataCleaner):
    """Clean the medata from Bruce dataset.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos
        """
        # Drop the "Unnamed: 0" column
        x.drop(columns=["Unnamed: 0"], inplace=True)

        # Add sample ID as index
        x.set_index("sample_id", inplace=True)

        # Replace Unknown values by nan
        x.replace("unknown", np.nan, inplace=True)

        # Remove columns with only nan values
        x.dropna(axis=1, how="all", inplace=True)

        return x


class BruceMIBIImagesCleaner(BaseDataCleaner):
    """Clean MIBI data file from Bruce.

    Add sample ID, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos, parsed dates
        """
        x.loc[:, "sample_id"] = x["path"].apply(lambda v: Path(v).parts[6])
        x.loc[:, "marker"] = x["path"].apply(lambda v: Path(v).stem)
        x = x.set_index("sample_id")
        return x


class BaseDataLoader:
    """Base class for loader."""

    def __init__(self, data_transformer=None, **load_data_kwargs):
        """Initialize base class.

        Parameters
        ----------
        data_transformer : optional
            An instance of DataCleaner or its subclass for data
            preprocessing, by default None
        **load_data_kwargs : dict
            Additional keyword arguments for the load_data method.
        """
        self.data_transformer = data_transformer
        self.load_data_kwargs = load_data_kwargs

    def load_data(self, data_path):
        """Load and transform the input data.

        Parameters
        ----------
        data_path : str
            Path to the data file to be read.

        Returns
        -------
        data : DataFrame
            Loaded and transformed data.
        """
        pass


class XLSXDataLoader(BaseDataLoader):
    """A class for loading Excel data.

    This transformer loads and preprocesses data using pandas.
    """

    def __init__(self, data_transformer=None, **load_data_kwargs):
        """Initialize the CSVDataLoader object.

        Parameters
        ----------
        data_transformer : optional
            An instance of DataCleaner for data
            preprocessing, default is None.
        load_data_kwargs : dict
            Additional keyword arguments for the load_data method.
        """
        # if data_transformer is None the data will be returned in the same
        # format as loaded

        super().__init__(data_transformer=data_transformer, **load_data_kwargs)

    def load_data(self, data_path):
        """Load and transform the input data.

        Parameters
        ----------
         data_path : str
            Path to the Excel file to be read.

        Returns
        -------
        X_csv : DataFrame
            Loaded and transformed data not including the target columns
        y : DataFrame
            Data corresponding to target columns with index as patient ids.
        """
        data = pd.read_excel(data_path, **self.load_data_kwargs)
        if self.data_transformer is not None:
            data = self.data_transformer.transform(data)
        return data


class CSVDataLoader(BaseDataLoader):
    """A class for loading CSV data.

    This transformer loads and preprocesses data using pandas.
    """

    def __init__(self, data_transformer=None, **load_data_kwargs):
        """Initialize the CSVDataLoader object.

        Parameters
        ----------
        dataset_name : str, optional
            Dataset name,
        data_transformer : optional
            An instance of DataCleaner or its subclass for data
            preprocessing, default is None.
        load_data_kwargs : dict
            Additional keyword arguments for the load_data method.
        """
        # if data_transformer is None the data will be returned in the same
        # format as loaded
        super().__init__(data_transformer=data_transformer, **load_data_kwargs)

    def load_data(
        self,
        data_path,
    ):
        """Transform the input data using pandas.

        Returns
        -------
        X_csv : DataFrame
            Loaded and transformed data not including the target columns
        y : DataFrame
            Data corresponding to target columns with index as patient ids.
        """
        data = pd.read_csv(filepath_or_buffer=data_path, **self.load_data_kwargs)
        if self.data_transformer is not None:
            data = self.data_transformer.transform(data)
        return data


class PathsLoader(BaseDataLoader):
    """Class to load a data frame with the path to the paths of files in a directory."""

    def __init__(
        self, pattern_files: str = "*", data_transformer=None, **load_data_kwargs
    ):
        """Initialize the PathsLoader object.

        Parameters
        ----------
        pattern_files : str
            Pattern to match the files in the directory.
        data_transformer : optional
            An instance of DataCleaner for data
            preprocessing, default is None.
        load_data_kwargs : dict
            Additional keyword arguments for the load_data method.
        """
        super().__init__(data_transformer=data_transformer, **load_data_kwargs)
        self.pattern_files = pattern_files

    def load_data(self, data_path):
        """Load and transform the input data.

        Parameters
        ----------
        data_path : str
            Path to the folder containing the files to be included.

        Returns
        -------
        data : DataFrame
            Loaded and transformed data.
        """
        paths = list(data_path.glob(self.pattern_files))
        data = pd.DataFrame({"path": paths})
        if self.data_transformer is not None:
            data = self.data_transformer.transform(data)
        return data


class SingleCellLoader:
    """Class to load the single cell anndata."""

    def load_data(self, path_data):
        """Load the single cell data.

        Parameters
        ----------
        path_data : str
            Path to the anndata object.

        Returns
        -------
        adata : anndata
            The anndata object.
        """
        adata = ad.read_h5ad(path_data)
        return adata


class VisiumLoader:
    """Class to load the visium anndata."""

    def load_samples_to_anndata(
        self,
        input_path,
        sample_list,
        resolution,
    ):
        """Load AnnData objects from a folder containing one folder per sample using scanpy.read_visium.

        Parameters
        ----------
        input_path (str): Path to the root folder containing sample subfolders.
                        Each subfolder should contain a Visium dataset with spatial data.
        sample_list (Optional[List[str]]): List of sample IDs to load. If None, all samples in the folder are loaded.
        resolution (str): Resolution of the spatial image to load, either "hires" or "lowres". Default is "hires".

        Returns
        -------
        Dict[str, sc.AnnData]: A dictionary where keys are sample IDs (folder names) and
                            values are the loaded AnnData objects.
        """
        if resolution not in ["hires", "lowres"]:
            raise ValueError("Resolution must be either 'hires' or 'lowres'.")

        anndata_dict = {}

        # Get the list of samples to process
        samples_to_load = (
            sample_list if sample_list is not None else os.listdir(input_path)
        )

        # Iterate over the specified samples
        for sample_id in samples_to_load:
            sample_path = os.path.join(input_path, sample_id, "outs")

            if os.path.isdir(sample_path):
                try:
                    # Load the Visium dataset using scanpy.read_visium
                    adata = sc.read_visium(sample_path)
                    adata.var_names_make_unique()
                    # Extract the spatial image and ensure the desired resolution is available
                    library_id = list(adata.uns["spatial"].keys())[0]
                    if resolution in adata.uns["spatial"][library_id]["images"]:
                        adata.uns["spatial"][library_id]["use_for_plotting"] = (
                            resolution
                        )
                    else:
                        print(
                            f"Resolution '{resolution}' not available for sample {sample_id}. Skipping."
                        )
                        continue

                    # Add the AnnData object to the dictionary
                    anndata_dict[sample_id] = adata

                except Exception as e:
                    print(f"Failed to load data for sample {sample_id}: {e}")
            else:
                print(f"Sample directory not found: {sample_path}")

        return anndata_dict


class DataFile(NamedTuple):
    """Define the structure of a datafile."""

    name: Optional[Union[str, Path]] = None
    loader: Optional[BaseDataLoader] = None

    def load(self):
        """Load datafile.

        Returns
        -------
        pd.DataFrame
            Loaded and cleaned datafile for a given source.

        Raises
        ------
        ValueError
            In case no loader has been specified.
        """
        if self.loader is not None:
            return self.loader.load_data(self.name)
        raise ValueError("Please specify a loader to load data.")


class DataSheet(NamedTuple):
    """Define the structure of an excel sheet."""

    name: str
    header: Optional[tuple[int]] = None


class DataSource(NamedTuple):
    """Define the structure of a data source (tabular modality)."""

    files: Union[DataFile, dict[str, DataFile]]


CENTERS: dict[str, DataCenter] = {}


def register_center(cls):
    """Define decorator to register instances of DataCenter."""

    @functools.wraps(cls)
    def wrapper_register(*args, **kwargs):
        obj = cls(*args, **kwargs)
        CENTERS[obj.name] = obj
        return obj

    return wrapper_register


@register_center
class DataCenter(NamedTuple):
    """Define the structure of a data center."""

    name: str
    sources: dict[str, DataSource]

    def load_tabular(self):
        """Load tabular data.

        Returns
        -------
        dict[str, pd.DataFrame]
            Dictionary of dataframes for each source that are tabular.
        """
        return {
            source_name: {
                file_name: (file_.load())
                for file_name, file_ in source.files.items()  # Only call items() if source.files is a dict
            }
            if isinstance(source.files, dict)
            else source.files.load()
            # Load only the sources that are clinical, bulk_rna
            for source_name, source in self.sources.items()
            if source_name in ["metadata", "clinical", "bulk_rna", "wes", "he", "mibi_images"]
        }

    def load_singlecell(self):
        """Load single-cell data.

        Returns
        -------
        dict[str, pd.DataFrame]
            Dictionary of dataframes for each source that are single-cell.
        """
        file_path = self.sources["sc_rna"].files["scRNA anndata"].name
        return SingleCellLoader().load_data(file_path)

    def load_visium(
        self, sample_list: Optional[List[str]] = None, resolution: str = "lowres"
    ):
        """Load Visium data.

        Parameters
        ----------
        sample_list (Optional[List[str]]): List of sample IDs to load. If None, all samples in the folder are loaded.
        resolution (str): Resolution of the spatial image to load, either "hires" or "lowres". Default is "hires".

        Returns
        -------
        Dict[str, sc.AnnData]: A dictionary where keys are sample IDs (folder names) and
                               values are the loaded AnnData objects.
        """
        spatial_files = self.sources["spatial"].files
        if isinstance(spatial_files, dict):
            visium_file = spatial_files["Visium anndata"]
        else:
            raise TypeError("Expected a dictionary for 'spatial' files")
        if isinstance(visium_file, DataFile):
            file_path = visium_file.name
        else:
            raise TypeError("Expected a DataFile for 'Visium anndata'")
        if sample_list is None:
            print("No sample list provided. Loading all samples.")
        print("Resolution of the spatial image to load: ", resolution)
        print(
            "You can change the resolution by setting the resolution parameter using the resolution argument."
        )
        print("Loading Visium data, this can take few minutes...")
        return VisiumLoader().load_samples_to_anndata(
            file_path, sample_list, resolution
        )


HOME_PATH = Path("/mnt/custom-file-systems/efs/fs-09913c1f7db79b6fd")
DATA_PATH = HOME_PATH / "data"
MOSAIC_PATH = DATA_PATH / "mosaic_dataset"

def generate_data_paths_mosaic(
    clinical_path: Path,
    bulkrna_path: Path,
    wes_path: Path,
    single_cell_path: Path,
    he_path: Path,
    he_features_path: Path,
    visium_path: Path,
):
    """Generate a dictionary of data paths."""
    return {
        "Processed clinical": clinical_path / "GBM_HK_sample_and_clinical_data.csv",
        "Key events clinical": clinical_path / "GBM_HK_multi_entry_events.csv",
        "Treatments": clinical_path / "GBM_HK_multi_entry_treatments.csv",
        "Data dictionary": clinical_path / "GBM_HK_data_dictionary.csv",
        "Bulk RNA raw counts": bulkrna_path / "counts/raw_counts.tsv",
        "Bulk RNA normalized counts": bulkrna_path / "counts/deseq2norm_counts.tsv",
        "Bulk RNA TPM counts": bulkrna_path / "counts/tpm_counts.tsv",
        "Bulk RNA fpkm counts": bulkrna_path / "counts/fpkm_counts.tsv",
        "WES CNV deletion": wes_path
        / "results/cohort_level/cnv/binary_deletion_status.csv",
        "WES CNV amplification": wes_path
        / "results/cohort_level/cnv/binary_duplication_status.csv",
        "WES CNV oncogenic": wes_path
        / "results/cohort_level/cnv/binary_oncogenic_alteration_status.csv",
        "WES mutations": wes_path
        / "results/cohort_level/snv_indel/binary_oncogenic_alteration_status.csv",
        "Single cell anndata": single_cell_path
        / "preprocessed/sc_merged_annotated.h5ad",
        "HE": he_path,
        "HE_features_H1": he_features_path,
        "Visium anndata": visium_path,
    }

DATA_PATH_FILES_MOSAIC = generate_data_paths_mosaic(
    MOSAIC_PATH / "Clinical",
    MOSAIC_PATH / "RNAseq",
    MOSAIC_PATH / "WES",
    MOSAIC_PATH / "single_cell",
    MOSAIC_PATH / "Visium" / "converted_he",
    DATA_PATH / "h1_bioptimus_features",  # H1 features
    MOSAIC_PATH / "Visium" / "spaceranger_count",
)


@frozen
class SourceNames:
    """Define the names of the sources.

    This class is used to avoid typos when accessing the sources of a center.
    """

    metadata: str = "metadata"
    clinical: str = "clinical"
    wes: str = "wes"
    bulk_rna: str = "bulk_rna"
    sc_rna: str = "sc_rna"
    spatial: str = "spatial"
    he: str = "he"
    mibi_images: str = "mibi_images"

    # force default values by forbidding arguments during initialization.
    def __init__(self):
        self.__attrs_init__()  # type: ignore # pylint: disable=no-member


SOURCES = SourceNames()

MosaicDataset = DataCenter(
    name="mosaic",
    sources={
        SOURCES.clinical: DataSource(
            files={
                "data dictionary": DataFile(
                    name=DATA_PATH_FILES_MOSAIC["Data dictionary"],
                    loader=CSVDataLoader(),
                ),
                "original clinical": DataFile(
                    name=DATA_PATH_FILES_MOSAIC["Processed clinical"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICLightClinicalDataCleaner(),
                    ),
                ),
                "processed gbm clinical": DataFile(
                    name=DATA_PATH_FILES_MOSAIC["Processed clinical"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICClinicalDataCleaner(),
                    ),
                ),
                "treatments": DataFile(
                    name=DATA_PATH_FILES_MOSAIC["Treatments"],
                    loader=CSVDataLoader(),
                ),
                "key events clinical": DataFile(
                    name=DATA_PATH_FILES_MOSAIC["Key events clinical"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICKeyEventsDataCleaner(),
                    ),
                ),
            }
        ),
        SOURCES.bulk_rna: DataSource(
            files={
                "raw counts": DataFile(
                    name=DATA_PATH_FILES_MOSAIC["Bulk RNA raw counts"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICBulkRNADataCleaner(), sep="\t"
                    ),
                ),
                "TPM counts": DataFile(
                    name=DATA_PATH_FILES_MOSAIC["Bulk RNA TPM counts"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICBulkRNADataCleaner(), sep="\t"
                    ),
                ),
                "normalized counts": DataFile(
                    name=DATA_PATH_FILES_MOSAIC["Bulk RNA normalized counts"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICBulkRNADataCleaner(), sep="\t"
                    ),
                ),
                "fpkm counts": DataFile(
                    name=DATA_PATH_FILES_MOSAIC["Bulk RNA fpkm counts"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICBulkRNADataCleaner(), sep="\t"
                    ),
                ),
            }
        ),
        SOURCES.spatial: DataSource(
            files={
                "Visium anndata": DataFile(
                    name=DATA_PATH_FILES_MOSAIC["Visium anndata"],
                    loader=None,
                ),
            }
        ),
        SOURCES.sc_rna: DataSource(
            files={
                "scRNA anndata": DataFile(
                    name=DATA_PATH_FILES_MOSAIC["Single cell anndata"],
                    loader=None,
                ),
            }
        ),
        SOURCES.wes: DataSource(
            files={
                "WES CNV deletion": DataFile(
                    name=DATA_PATH_FILES_MOSAIC["WES CNV deletion"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICWESDataCleaner(),
                    ),
                ),
                "WES CNV amplification": DataFile(
                    name=DATA_PATH_FILES_MOSAIC["WES CNV amplification"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICWESDataCleaner(),
                    ),
                ),
                "WES CNV oncogenic": DataFile(
                    name=DATA_PATH_FILES_MOSAIC["WES CNV oncogenic"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICWESDataCleaner(),
                    ),
                ),
                "WES mutations": DataFile(
                    name=DATA_PATH_FILES_MOSAIC["WES mutations"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICWESDataCleaner(),
                    ),
                ),
            }
        ),
        SOURCES.he: DataSource(
            files={
                "HE files": DataFile(
                    name=DATA_PATH_FILES_MOSAIC["HE"],
                    loader=PathsLoader(
                        pattern_files="*.tif",
                        data_transformer=MOSAICHEDataCleaner(),
                    ),
                ),
                "H1 features": DataFile(
                    name=DATA_PATH_FILES_MOSAIC["HE_features_H1"],
                    loader=PathsLoader(
                        pattern_files="*_eHnE.zarr",
                        data_transformer=MOSAICHEFeaturesDataCleaner(),
                    ),
                ),
            },
        ),
    },
)
