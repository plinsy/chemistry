from chemistry import generer_tous_squelettes, calculer_insaturation

# Test C6H6
c = 6
h = 6
dou = calculer_insaturation(c, h)

print(f"C{c}H{h}")
print(f"Degré d'insaturation: {dou}")
print(f"\nGénération des isomères...")

resultats = generer_tous_squelettes(c, h)

print(f"\nNombre d'isomères trouvés: {len(resultats)}")
print(f"Attendu: 217")
print(f"Différence: {217 - len(resultats)}")

if len(resultats) < 30:
    print("\nListe des SMILES:")
    for idx, smiles in enumerate(sorted(resultats), 1):
        print(f"  {idx}. {smiles}")
