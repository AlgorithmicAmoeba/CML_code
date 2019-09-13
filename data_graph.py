import pandas
import matplotlib.pyplot as plt
import itertools

conc7 = pandas.read_csv("data/run_7_conc.csv")
glucose7 = pandas.read_csv("data/run_7_glucose.csv")

conc9 = pandas.read_csv("data/run_9_conc.csv")
glucose9 = pandas.read_csv("data/run_9_glucose.csv")

plt.figure(figsize=(20, 40))
for i, (conc, glucose) in enumerate([(conc7, glucose7), (conc9, glucose9)]):
    plt.subplot(2, 4, 4*i + 1)
    plt.plot(glucose['Time'], glucose['Glucose dosing (g/h)'])
    plt.title("Glucose feed (g/h)")

    for j, name in enumerate(["Glucose", "Fumaric", "Ethanol"]):
        plt.subplot(2, 4, 4*i + j + 2)
        plt.plot(conc['Time'], conc[name], '.')
        plt.title(f"Concentration of {name} (g/L)")

plt.show()

plt.figure(figsize=(20, 40))
for i, (conc, glucose) in enumerate([(conc7, glucose7), (conc9, glucose9)]):
    combos = itertools.combinations(['Glucose', 'Fumaric', 'Ethanol'], 2)
    for j, combo in enumerate(combos):
        a, b = combo
        plt.subplot(2, 3, 3*i + j + 1)
        plt.plot(conc[a], conc[b], '.')
        plt.xlabel(a)
        plt.ylabel(b)

plt.show()

plt.figure(figsize=(20, 40))
for i, (conc, glucose) in enumerate([(conc7, glucose7), (conc9, glucose9)]):
    combos = itertools.permutations(['Glucose', 'Fumaric', 'Ethanol'], 3)
    for j, combo in enumerate(combos):
        a, b, c = combo
        plt.subplot(2, 6, 6*i + j + 1)
        plt.plot(conc[a], conc[b]/conc[c], '.')
        plt.xlabel(a)
        plt.ylabel(b + ' / ' + c)

plt.show()
