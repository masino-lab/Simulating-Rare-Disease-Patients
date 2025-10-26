import random, csv

#set seed to be reproduceable
random.seed(42)

#define the true phenotypes associated with each of the 5 diseases in a list
TRUE_PHENOTYPES = {
    #Multiple epiphyseal dysplasia, Al-Gazali type
    "D001": ["HP:0000256", "HP:0000272", "HP:0000316", "HP:0000369", "HP:0000470", "HP:0000767",
        "HP:0001274", "HP:0001373", "HP:0001513", "HP:0002007", "HP:0002758", "HP:0002857", 
        "HP:0005930", "HP:0006101", "HP:0012444", "HP:0030084"],
    #Beta-mannosidosis
    "D002": ["HP:0000365", "HP:0001249", "HP:0001250", "HP:0001999", "HP:0002205", "HP:0005247"],
    #Glycogen storage disease due to muscle phosphofructokinase deficiency
    "D003": ["HP:0001324", "HP:0001903", "HP:0002149", "HP:0002486"],
    #Isolated osteopoikilosis
    "D004": ["HP:0000086", "HP:0000252", "HP:0001482", "HP:0002652", "HP:0004322", "HP:0005789"],
    #Hyperekplexia-epilepsy syndrome
    "D005": ["HP:0001276", "HP:0002267", "HP:0002376", "HP:0002384", "HP:0007333", 
        "HP:0010818", "HP:0012018", "HP:0200134", "HP:0000243"]
}

#create a noisy array (20 noisy HPO)
NOISY_HPOS = ["HP:0001990", "HP:0002024", "HP:0002167", "HP:0002204", "HP:0002240", 
    "HP:0002360", "HP:0002650", "HP:0002684", "HP:0002750", "HP:0002997", "HP:0003103", 
    "HP:0003196", "HP:0003468", "HP:0004337", "HP:0004568", "HP:0008430", "HP:0008551", 
    "HP:0011276", "HP:0012068", "HP:0012471"]

rows = []
for pid in range (1,101): # make 100 mock patients

    #assign a random disease from the list above
    disease = random.choice(list(TRUE_PHENOTYPES.keys()))
    true_terms = TRUE_PHENOTYPES[disease]

    # start with all of the true terms
    observed = set(true_terms)

    #randomly drop some of the true terms (false negative)
    n_missing = random.randint(0,2)
    if n_missing:
        dropped = set(random.sample(true_terms, n_missing))
    else:
        dropped = set()
    
    for n in dropped:
        observed.discard(n)

    # randomly add noisy phenotypes
    n_noisy = random.randint(0,3)
    added_noise = random.sample(NOISY_HPOS, n_noisy)

    # write the observed phenotypes in rows
    for hpo in observed:
        rows.append([pid, disease, hpo, 1])
    for hpo in added_noise:
        rows.append([pid, disease, hpo, 0])

# write this out
with open("mock_phenotypes.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["patient_id", "disease_id", "hpo_id", "is_true"])
    w.writerows(rows)

print("Wrote mock_phenotypes.csv with", len(rows), "rows.")

