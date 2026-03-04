from flask import Flask, render_template, request, jsonify
from flask_htmx import HTMX
from chemistry import generer_tous_squelettes, parse_formule
from rdkit import Chem
from rdkit.Chem import AllChem
from collections import OrderedDict

app = Flask(__name__)
htmx = HTMX(app)

PAGE_SIZE = 20
MAX_CACHE_ENTRIES = 50
MAX_VIEWED_FORMULAS = 10

FORMULA_CACHE = OrderedDict()
VIEWED_FORMULAS = OrderedDict()


def molecule_to_view_data(mol):
    conf = mol.GetConformer()
    atoms_data = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        atoms_data.append(
            {
                "element": atom.GetSymbol(),
                "x": pos.x,
                "y": pos.y,
                "z": pos.z,
            }
        )

    bonds_data = []
    for bond in mol.GetBonds():
        bonds_data.append(
            {
                "start": bond.GetBeginAtomIdx(),
                "end": bond.GetEndAtomIdx(),
                "type": float(bond.GetBondTypeAsDouble()),
            }
        )

    return atoms_data, bonds_data


def _remember_viewed_formula(formula):
    VIEWED_FORMULAS[formula] = True
    VIEWED_FORMULAS.move_to_end(formula)

    while len(VIEWED_FORMULAS) > MAX_VIEWED_FORMULAS:
        VIEWED_FORMULAS.popitem(last=False)


def _build_isomer_payload(smiles_list):
    isomers = []
    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)  # type: ignore

            atoms_data, bonds_data = molecule_to_view_data(mol)

            isomers.append({"smiles": smiles, "atoms": atoms_data, "bonds": bonds_data})
        except Exception:
            try:
                fallback_mol = Chem.MolFromSmiles(smiles)
                fallback_mol = Chem.AddHs(fallback_mol)
                AllChem.Compute2DCoords(fallback_mol)
                atoms_data, bonds_data = molecule_to_view_data(fallback_mol)
                isomers.append({"smiles": smiles, "atoms": atoms_data, "bonds": bonds_data})
            except Exception:
                isomers.append({"smiles": smiles, "atoms": [], "bonds": []})

    return isomers


@app.route("/")
def home():
    if htmx:
        return render_template("partials/thing.html")
    return render_template("pages/index.html")


@app.route("/api/calculate", methods=["POST"])
def calculate_isomers():
    """API endpoint to calculate molecular isomers"""
    try:
        # Get formula from form data
        formula = request.form.get("formula", "").strip().upper()

        if not formula:
            return (
                render_template(
                    "partials/error.html", error="Please enter a chemical formula"
                ),
                400,
            )

        # Parse the formula
        atoms = parse_formule(formula)

        # Extract carbon and hydrogen counts
        c = atoms.get("C", 0)
        h = atoms.get("H", 0)

        if c == 0:
            return (
                render_template(
                    "partials/error.html", error="Formula must contain carbon atoms"
                ),
                400,
            )

        page = request.form.get("page", "1").strip()
        try:
            current_page = max(1, int(page))
        except ValueError:
            current_page = 1

        is_cached = formula in FORMULA_CACHE
        if is_cached:
            FORMULA_CACHE.move_to_end(formula)
            isomers = FORMULA_CACHE[formula]
        else:
            smiles_list = generer_tous_squelettes(c, h)
            isomers = _build_isomer_payload(smiles_list)
            FORMULA_CACHE[formula] = isomers
            FORMULA_CACHE.move_to_end(formula)

            while len(FORMULA_CACHE) > MAX_CACHE_ENTRIES:
                FORMULA_CACHE.popitem(last=False)

        _remember_viewed_formula(formula)

        total_isomers = len(isomers)
        total_pages = max(1, (total_isomers + PAGE_SIZE - 1) // PAGE_SIZE)
        current_page = min(current_page, total_pages)
        start_idx = (current_page - 1) * PAGE_SIZE
        end_idx = min(start_idx + PAGE_SIZE, total_isomers)
        paginated_isomers = isomers[start_idx:end_idx]

        page_ranges = []
        if total_isomers > 0:
            for page_number in range(1, total_pages + 1):
                range_start = (page_number - 1) * PAGE_SIZE + 1
                range_end = min(page_number * PAGE_SIZE, total_isomers)
                page_ranges.append(
                    {
                        "number": page_number,
                        "start": range_start,
                        "end": range_end,
                    }
                )

        return render_template(
            "partials/molecules.html",
            isomers=paginated_isomers,
            formula=formula,
            total_isomers=total_isomers,
            current_page=current_page,
            total_pages=total_pages,
            start_index=start_idx + 1 if total_isomers > 0 else 0,
            end_index=end_idx,
            page_ranges=page_ranges,
            viewed_formulas=list(reversed(VIEWED_FORMULAS.keys())),
            is_cached=is_cached,
        )

    except Exception as e:
        return (
            render_template(
                "partials/error.html", error=f"Error calculating isomers: {str(e)}"
            ),
            500,
        )
