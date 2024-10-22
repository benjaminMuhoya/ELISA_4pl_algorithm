#!/usr/bin/env python3
import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from collections import defaultdict
import seaborn as sns
import pandas as pd

def four_parameter_logistic(x, A, B, C, D):
    """4PL logistic equation."""
    return ((A-D)/(1.0+(np.sign((x/C))*np.abs(x/C)**B))) + D

def four_parameter_logistic_inverse(g, A, B, C, D):
    return C*np.sign((A-D)/(g-D)-1)*np.abs((A-D)/(g-D)-1)**(1/B)

# Read in data
with open('last_pilot.csv', 'r', encoding='utf-8-sig') as f:
    reader = csv.reader(f)
    data = list(reader)

# Separate columns
sample_names = [i[0] for i in data]
sample_x_values = [i[1] for i in data]
standard_x_values1 = [i[2] for i in data]
standard_x_values2 = [i[3] for i in data]
standard_y_values = [i[4] for i in data]
blank_x_values = [i[5] for i in data]

# Remove blanks
sample_x_values = [x for x in sample_x_values if x]
standard_x_values1 = [x for x in standard_x_values1 if x]
standard_x_values2 = [x for x in standard_x_values2 if x]
standard_y_values = [x for x in standard_y_values if x]
blank_x_values = [x for x in blank_x_values if x]

# Convert to floats
sample_x_values = [float(i) for i in sample_x_values]
standard_x_values1 = [float(i) for i in standard_x_values1]
standard_x_values2 = [float(i) for i in standard_x_values2]
standard_y_values = [float(i) for i in standard_y_values]
blank_x_values = [float(i) for i in blank_x_values]

print("Raw Sample Absorbance Values: " + str(sample_x_values))
print("Standard Replicate1 Raw Absorbance Values: " + str(standard_x_values1))
print("Standard Replicate2 Raw Absorbance Values: " + str(standard_x_values2))
print("Standard Protein Concentration Values (ng/mL): " + str(standard_y_values))
print("Blank Well Raw Absorbance Values: " + str(blank_x_values))

# Get averages for R
R = []
for i in range(len(standard_x_values1)):
    R.append((standard_x_values1[i]+standard_x_values2[i])/2)

print("Average Standard Absorbance Values: " + str(R))

# Get average for blank x values
blank_x_average = sum(blank_x_values)/len(blank_x_values)

print("Blank Correction Factor: " + str(blank_x_average))

# Get W
W = []
for i in range(len(R)):
    W.append(R[i] - blank_x_average)

print("Blank-corrected Standard absorbance values: " + str(W))

# Get Z
Z = []
for i in range(len(sample_x_values)):
    Z.append(sample_x_values[i] - blank_x_average)

print("Blank-corrected Sample absorbance values: " + str(Z))

# Fit curve
popt, pcov = curve_fit(four_parameter_logistic, standard_y_values, W, maxfev=20000)
A, B, C, D = popt

print("Parameter Values (A, B, C, D):" + str(popt))

# Get results in ng/mL
results_ng_ml = four_parameter_logistic_inverse(Z, A, B, C, D)

# Calculate actual concentration in original sample after 1500-fold dilution (convert to µg/mL)
dilution_factor = 1500
results_ug_ml = [(conc * dilution_factor) / 1000 for conc in results_ng_ml]  # Convert from ng/mL to µg/mL

print("Protein concentrations of samples in ng/mL: " + str(results_ng_ml))
print("Protein concentrations of samples in µg/mL (after 1500-fold dilution): " + str(results_ug_ml))

# Get R-squared
residuals = W - four_parameter_logistic(standard_y_values, *popt)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((W-np.mean(W))**2)
r_squared = 1 - (ss_res / ss_tot)
print("R-Squared Value: " + str(r_squared))

# Output results
with open('last_pilot_OUTPUT.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['This 4-Parameter Logistic Curve fit output was generated by Genesis Lung in '
                     'tandem with ChatGPT (03/02/2023 Version).'])
    writer.writerow([])
    writer.writerow(['4PL Logistic Equation', '4PL Parameter A', '4PL Parameter B', '4PL Parameter C',
                     '4PL Parameter D'])
    writer.writerow(['D+(A-D)/(1.0+(x/C)^B)', A, B, C, D])
    writer.writerow([])
    writer.writerow(['R-Squared Value', 'Blank Correction Factor'])
    writer.writerow([r_squared, blank_x_average])
    writer.writerow([])
    writer.writerow(['Blank Well Absorbance Values'])
    for i in range(len(blank_x_values)):
        writer.writerow([blank_x_values[i]])
    writer.writerow([])
    writer.writerow(['Raw Replicate 1 Standard Absorbance Values', 'Raw Replicate 2 Standard Absorbance Values',
                     'Average Blank-Corrected Standard Absorbance', 'Standard Protein Conc. (pg/mL)'])
    for i in range(len(W)):
        writer.writerow([standard_x_values1[i], standard_x_values2[i], W[i], standard_y_values[i]])
    writer.writerow([])
    writer.writerow(['Sample ID', 'Raw Sample Absorbance Values', 'Blank-Corrected Sample Absorbance',
                     'Sample Concentration (ng/mL)', 'Adjusted Concentration (µg/mL after dilution)'])
    for i in range(len(sample_x_values)):
        writer.writerow([sample_names[i], sample_x_values[i], Z[i], results_ng_ml[i], results_ug_ml[i]])

# Plot Standard Curve
xaxis = np.linspace(0, 160)
yaxis = ((A-D)/(1.0+(np.sign((xaxis/C))*np.abs(xaxis/C)**B))) + D
baxis = []
for i in range(len(blank_x_values)):
    baxis.append(0)
plt.figure()
plt.title('LBP Concentration VS Absorbance')
plt.plot(xaxis, yaxis, 'r-', label='4-Parameter Logistic Curve Fit (r2={:.3f})'.format(r_squared))
plt.plot(results_ng_ml, Z, 'go', label='Experimental Samples')
plt.plot(standard_y_values, W, 'ro', label='Averaged Standards')
plt.plot(standard_y_values, standard_x_values1, 'yo', label='Standard Replicate 1')
plt.plot(standard_y_values, standard_x_values2, 'bo', label='Standard Replicate 2')
plt.plot(baxis, blank_x_values, 'ko', label='Blank Wells')
plt.legend()
plt.xlabel('Absolute Conc. (ng/mL)')
plt.ylabel('Blank-Corrected Absorbance')
plt.show()

# --- New Plot for Sample Concentrations with CV ---
# Group results based on identical sample names
sample_groups = defaultdict(list)
for i, name in enumerate(sample_names):
    sample_groups[name].append(results_ug_ml[i])

# Filter out samples that don't have duplicates
duplicate_samples = {k: v for k, v in sample_groups.items() if len(v) > 1}
num_duplicates = len(duplicate_samples)

# Print a message with the number of duplicate sample names found
print(f"Number of duplicate sample names found: {num_duplicates}")

if num_duplicates > 0:
    plt.figure(figsize=(8, 6))
    
    # Plot the concentrations and calculate CV for each duplicate sample group
    for sample_name, concentrations in duplicate_samples.items():
        plt.scatter([sample_name] * len(concentrations), concentrations, label=sample_name, color='blue')

        # Calculate the mean and standard deviation for CV
        mean_conc = np.mean(concentrations)
        std_conc = np.std(concentrations)
        cv = (std_conc / mean_conc) * 100 if mean_conc != 0 else 0

        # Display CV on the plot above each sample group
        plt.text(sample_name, max(concentrations) + 0.5, f'CV: {cv:.2f}%', ha='center', color='black', fontsize=12, rotation=45)

    plt.title('Sample Concentrations with CV (in µg/mL)')
    plt.xlabel('Sample Name')
    plt.ylabel('Concentration (µg/mL)')
    plt.grid(True)
    plt.xticks(rotation=45, ha='right')
    plt.show()
else:
    print("No duplicate sample names found.")

# --- Third Plot: Distribution of Concentrations (Boxplot + Violin Plot) ---
# Convert results to a DataFrame for easier plotting
df = pd.DataFrame({
    'Sample_Name': sample_names,
    'Concentration_ug_ml': results_ug_ml
})

# Plot the distribution with a boxplot overlaid with a violin plot
plt.figure(figsize=(8, 6))
sns.violinplot(y='Concentration_ug_ml', data=df, inner=None, color='lightblue')
sns.boxplot(y='Concentration_ug_ml', data=df, whis=1.5, width=0.2, color='orange')
plt.title('Distribution of Adjusted Concentration (µg/mL)')
plt.grid(True)
plt.show()

