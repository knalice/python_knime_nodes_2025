import knime.extension as knext
from rdkit import Chem
import logging
LOGGER = logging.getLogger(__name__)

category = knext.category(
    '/community/rdkit',
    'python_nodes',
    'New RDKit Nodes',
    'RDKit nodes which are written in Python (still exploratory)',
    icon='category_rdkit.png')


_ctabTypes = None
_smilesTypes = None
_rdkitTypes = None
def _setup_knime_types():
    global _ctabTypes, _smilesTypes, _rdkitTypes
    if _ctabTypes is None:
        import knime.types.chemistry as cet  # To work with and compare against chemical data types like SMILES,...
        _ctabTypes = (
            knext.logical(cet.SdfValue),
            knext.logical(cet.SdfAdapterValue),
            knext.logical(cet.MolAdapterValue),
        )
        _smilesTypes = (
            knext.logical(cet.SmilesValue),
            knext.logical(cet.SmilesAdapterValue),
        )
        _rdkitTypes = (knext.logical(Chem.rdchem.Mol), )


def column_is_convertible_to_mol(column: knext.Column):
    _setup_knime_types()
    c_type = column.ktype
    allowedTypes = _smilesTypes + _ctabTypes + _rdkitTypes

    return c_type in allowedTypes


def convert_column_to_rdkit_mol(df,
                                molecule_column_type,
                                molecule_column_param,
                                sanitizeOnParse=True):
    _setup_knime_types()
    if molecule_column_type in _rdkitTypes:
        LOGGER.warning("rdkit mols detected")
        mols = df[molecule_column_param]
    elif molecule_column_type in _smilesTypes:
        mols = [
            Chem.MolFromSmiles(smi, sanitize=sanitizeOnParse)
            for smi in df[molecule_column_param]
        ]
    elif molecule_column_type in _ctabTypes:
        mols = [
            Chem.MolFromMolBlock(mb, sanitize=sanitizeOnParse, removeHs=sanitizeOnParse)
            for mb in df[molecule_column_param]
        ]
    else:
        raise ValueError('unrecognized molecule column type')
    return mols
