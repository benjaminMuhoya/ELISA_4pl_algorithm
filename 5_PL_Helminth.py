#!/usr/bin/env python3
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
from collections import defaultdict

# Define the quadratic model
def quadratic_model(x, a, b, c):
    return a * x**2 + b * x + c

# Load the data
file_path = 'Helminth_dilutions.csv'  # Ensure this file is in the same directory or provide the full path
with open(file_path, 'r', encoding='utf-8-sig') as f:
    reader = csv.reader(f)
    data = list(reader)

# Separate columns
sample_names = [i[0] for i in data]
sample_x_values = [i[1] for i in data]
standard_x_values1 = [i[2] for i in data]
standard_x_values2 = [i[3] for i in data]
standard_y_values = [i[4] for i in data]
blank_x_values = [i[5] for i in data]

# Remove blanks and convert to floats
sample_x_values = [float(x) for x in sample_x_values if x]
standard_x_values1 = [float(x) for x in standard_x_values1 if x]
standard_x_values2 = [float(x) for x in standard_x_values2 if x]
standard_y_values = [float(y) for y in standard_y_values if y]
blank_x_values = [float(x) for x in blank_x_values if x]

# Calculate R (average standard absorbance) and W (blank-corrected standard absorbance)
R = [(standard_x_values1[i] + standard_x_values2[i]) / 2 for i in range(len(standard_x_values1))]
blank_x_average = np.mean(blank_x_values)
W = [r - blank_x_average for r in R]
Z = [x - blank_x_average for x in sample_x_values]

# Fit the quadratic model to standard data
standard_conc = np.array(standard_y_values)
standard_abs = np.array(W)
quadratic_params, _ = curve_fit(quadratic_model, standard_conc, standard_abs)

# Define parameters
a, b, c = quadratic_params
print("Quadratic Model Parameters (a, b, c):", quadratic_params)

# Invert the quadratic model to infer concentrations from absorbance values
def quadratic_inverse(y, a, b, c):
    A = a
    B = b
    C = c - y
    discriminant = B**2 - 4*A*C
    if discriminant < 0:
        return np.nan  # No real solution
    else:
        return (-B + np.sqrt(discriminant)) / (2*A)  # Only positive root

# Calculate concentrations for each sample using the inverse function
results_ng_ml = [quadratic_inverse(z, a, b, c) for z in Z]
dilution_factor = 5000
results_ug_ml = [(conc * dilution_factor) / 1000 for conc in results_ng_ml]  # Convert from ng/mL to µg/mL
results_mg_dl = [conc * dilution_factor / 10000 for conc in results_ng_ml]  # Convert ng/mL to mg/dL

# Output results to CSV
output_file = 'Helminth_dilutions_OUTPUT_QUADRATIC_with_mg_dl.csv'
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Quadratic Model Fit Results'])
    writer.writerow([])
    writer.writerow(['Quadratic Model Equation', 'Parameter a', 'Parameter b', 'Parameter c'])
    writer.writerow(['ax^2 + bx + c', a, b, c])
    writer.writerow([])
    writer.writerow(['Blank Correction Factor'])
    writer.writerow([blank_x_average])
    writer.writerow([])
    writer.writerow(['Sample ID', 'Raw Sample Absorbance', 'Blank-Corrected Absorbance',
                     'Concentration (ng/mL)', 'Adjusted Concentration (µg/mL after dilution)',
                     'Adjusted Concentration (mg/dL after dilution)'])
    for i in range(len(sample_x_values)):
        writer.writerow([sample_names[i], sample_x_values[i], Z[i], results_ng_ml[i], results_ug_ml[i], results_mg_dl[i]])

# Plot standard curve with quadratic model and sample scatter in green
xaxis = np.linspace(0, max(standard_conc)*1.2, 100)
yaxis = quadratic_model(xaxis, *quadratic_params)
plt.figure()
plt.plot(xaxis, yaxis, label="Quadratic Fit", color='blue')
plt.scatter(standard_conc, W, color='red', label="Standard Data")
plt.scatter(results_ng_ml, Z, color='green', label="Sample Data")
plt.title("Standard Curve with Quadratic Fit")
plt.xlabel("Protein Concentration (mg/dL)")
plt.ylabel("Absorbance (Blank Corrected)")
plt.legend()
plt.grid(True)
plt.show()

# Group results by sample name and calculate CV for duplicates
sample_groups = defaultdict(list)
for i, name in enumerate(sample_names):
    sample_groups[name].append(results_mg_dl[i])

# Filter for duplicate samples and calculate CV
duplicate_samples = {k: v for k, v in sample_groups.items() if len(v) > 1}
if duplicate_samples:
    plt.figure(figsize=(8, 6))
    for sample_name, concentrations in duplicate_samples.items():
        plt.scatter([sample_name] * len(concentrations), concentrations, label=sample_name, color='blue')
        mean_conc = np.mean(concentrations)
        std_conc = np.std(concentrations)
        cv = (std_conc / mean_conc) * 100 if mean_conc != 0 else 0
        plt.text(sample_name, max(concentrations) + 0.5, f'CV: {cv:.2f}%', ha='center', color='black', fontsize=10)

    plt.title('Sample Concentrations with CV (in mg/dL)')
    plt.xlabel('Sample Name')
    plt.ylabel('Concentration (mg/dL)')
    plt.grid(True)
    plt.xticks(rotation=45, ha='right')
    plt.show()
else:
    print("No duplicate sample names found.")

# Violin + Boxplot for Distribution of Adjusted Concentrations in mg/dL
df = pd.DataFrame({'Sample_Name': sample_names, 'Concentration_mg_dl': results_mg_dl})
plt.figure(figsize=(8, 6))
sns.violinplot(y='Concentration_mg_dl', data=df, inner=None, color='lightblue')
sns.boxplot(y='Concentration_mg_dl', data=df, whis=1.5, width=0.2, color='orange')
plt.title('Distribution of Adjusted Concentration (mg/dL)')
plt.grid(True)
plt.show()

