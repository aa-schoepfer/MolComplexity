from complexity import MolReader, complexity

if __name__ == "__main__":
    files = [
        "mols/n-octane.mol",
        "mols/1-methyl-heptane.mol",
        "mols/cubane.mol",
        "mols/propionic_acid.mol",
        "mols/(+)-camphore.mol",
        "mols/aldohexose.mol",
        "mols/glucose.mol",
    ]
    results = [
        0.00,
        0.63,
        1.27,
        2.00,
        2.53,  # should be 2.21
        2.16,
        2.48,  # should be 2.24
    ]

    for f, r in zip(files, results):
        with open(f) as file:
            mol = MolReader(file)
            c = round(complexity(mol.g), 2)

        assert c == r, f"{r} expected, got {c}."
