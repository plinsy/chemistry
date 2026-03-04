from flask import Flask, render_template, request, jsonify
from flask_htmx import HTMX
from chemistry import generer_tous_squelettes, parse_formule
from rdkit import Chem
from rdkit.Chem import AllChem

app = Flask(__name__)
htmx = HTMX(app)


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

        # Generate all isomers
        smiles_list = generer_tous_squelettes(c, h)

        # Prepare isomers data with 3D coordinates including hydrogens
        isomers = []
        for smiles in smiles_list:
            try:
                # Convert SMILES to molecule and add explicit hydrogens
                mol = Chem.MolFromSmiles(smiles)
                mol = Chem.AddHs(mol)

                # Generate 3D coordinates
                AllChem.EmbedMolecule(mol, randomSeed=42)  # type: ignore

                atoms_data, bonds_data = molecule_to_view_data(mol)

                isomers.append(
                    {"smiles": smiles, "atoms": atoms_data, "bonds": bonds_data}
                )
            except Exception as e:
                try:
                    # Fallback: still provide topology with 2D coordinates
                    fallback_mol = Chem.MolFromSmiles(smiles)
                    fallback_mol = Chem.AddHs(fallback_mol)
                    AllChem.Compute2DCoords(fallback_mol)
                    atoms_data, bonds_data = molecule_to_view_data(fallback_mol)
                    isomers.append(
                        {"smiles": smiles, "atoms": atoms_data, "bonds": bonds_data}
                    )
                except Exception:
                    # Final fallback: include SMILES only
                    isomers.append({"smiles": smiles, "atoms": [], "bonds": []})

        # Return the molecules partial with isomers data
        return render_template(
            "partials/molecules.html", isomers=isomers, formula=formula
        )

    except Exception as e:
        return (
            render_template(
                "partials/error.html", error=f"Error calculating isomers: {str(e)}"
            ),
            500,
        )
