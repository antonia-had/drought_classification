import pandas as pd
from scipy.stats import genextreme as gev
from utils import *


# load paleo data at Cisco
Paleo = pd.read_csv('./data/Reconstruction/Cisco_Recon_v_Observed_v_Stateline.csv')

# re-scale Cisco data to estimate data at CO-UT state line
factor = np.nanmean(Paleo['ObservedNaturalStateline']/Paleo['ObservedNaturalCisco'])
Paleo['ScaledNaturalCisco'] = Paleo['ObservedNaturalCisco']*factor*1233.4818/1000000
Paleo['ScaledReconCisco'] = Paleo['ReconCisco']*factor*1233.4818/1000000

historic_flows = np.load('./data/historic_flows.npy')
annual_historic_flows = np.sum(historic_flows, axis=1)*1233.4818/1000000

mean_flow = np.mean(annual_historic_flows)
st_deviation = np.std(annual_historic_flows)
drought_threshold = mean_flow-0.5*st_deviation

# Droughts under historical observations and under paleo reconstruction
drought_identification_plots(Paleo['ScaledReconCisco'][:429], mean_flow, drought_threshold)
drought_identification_plots(annual_historic_flows, mean_flow, drought_threshold)

# Repeat for synthetically generated flows (internal variability of history)
stationary_flows = np.load('./data/stationarysynthetic_flows.npy')*1233.4818/1000000
annual_stationary_flows = np.sum(stationary_flows, axis=2)
annual_stationary_flows_flat = annual_stationary_flows.flatten()
drought_counts(annual_stationary_flows_flat, drought_threshold, decadal_window=11, multidecadal_window=35)
drought_identification_plots(annual_stationary_flows_flat, mean_flow, drought_threshold)

# Repeat for synthetically generated flows (internal variability of paleo)
paleo_flows = np.load('./data/Paleo_SOWs_flows.npy')
paleo_flows_annual = np.sum(paleo_flows, axis=2)*1233.4818/1000000
# Flatten every 10 realizations
paleo_flows_annual_flat = paleo_flows_annual[0*10:0*10+10].flatten()
for i in range(1, 366):
    paleo_flows_annual_flat = np.vstack((paleo_flows_annual_flat, paleo_flows_annual[i*10:i*10+10].flatten()))

multi_numbers = []
decadal_numbers = []
for j in range(len(paleo_flows_annual_flat[:,0])):
    [multi, decadal] = drought_counts(paleo_flows_annual_flat[j,:], drought_threshold, decadal_window=11, multidecadal_window=35)
    multi_numbers.append(multi)
    decadal_numbers.append(decadal)



# Repeat for all encompassing sample
generated_flows_wider = np.load('./data/LHsamples_wider_100_AnnQonly_flows.npy')
all_annual_experiment_flows = np.sum(generated_flows_wider, axis=2)*1233.4818/1000000
# Flatten every 10 realizations
all_annual_experiment_flat = all_annual_experiment_flows[0*10:0*10+10].flatten()
for i in range(1, 100):
    all_annual_experiment_flat = np.vstack((all_annual_experiment_flat, all_annual_experiment_flows[i*10:i*10+10].flatten()))

multi_numbers = []
decadal_numbers = []
for j in range(len(all_annual_experiment_flat[:,0])):
    [multi, decadal] = drought_counts(all_annual_experiment_flat[j,:], drought_threshold, decadal_window=11, multidecadal_window=35)
    multi_numbers.append(multi)
    decadal_numbers.append(decadal)

# Summarize impacts in bar chart
# Order is: History, Paleo, Synthetic History, Synthetic Paleo, All-encompassing
decadal_means, multi_means = (21, 12, 46, 26, 57), (0, 0, 20, 5, 51)
decadal_ranges = [[0, 0, 0, 21, 57],
                  [0, 0, 0, 31, 43]]
multi_ranges = [[0, 0, 0, 5, 51],
                [0, 0, 0, 37, 49]]

ind = np.arange(len(decadal_means))  # the x locations for the groups
width = 0.35  # the width of the bars

fig, ax = plt.subplots(figsize=(9,6))
rects1 = ax.bar(ind - width/2, decadal_means, width, yerr=decadal_ranges,
                label='Decadal')
rects2 = ax.bar(ind + width/2, multi_means, width, yerr=multi_ranges,
                label='Multidecadal')

ax.set_ylabel('# of years in each\n drought per century', fontsize=16)
ax.set_xticks(ind)
ax.set_xticklabels(('History', 'Paleo', 'Synthetic\nHistory', 'Synthetic\nPaleo', 'All-\nencompassing'))
ax.tick_params(axis='both', labelsize=12)
ax.legend(fontsize=14, loc='upper left')
plt.show()

zero_droughts = [i for i, e in enumerate(decadal_numbers) if e == 0]
probability_decadal = 100*(1-len(zero_droughts)/len(decadal_numbers))

zero_multi = [i for i, e in enumerate(multi_numbers) if e == 0]
probability_multi = 100*(1-len(zero_multi)/len(decadal_numbers))

realizations_of_1000 = [sum(decadal_numbers[i:i+10]) for i in range(0, len(decadal_numbers), 10)]

hist, bins = np.histogram(realizations_of_1000, bins=9, density=False)
# convert histogram values to percentages
hist = np.around(hist/len(realizations_of_1000)*100, decimals=0)
fig, (ax1, ax2) = plt.subplots(2, figsize=(6,9))

ax1.bar(bins[:-1], hist, width=np.diff(bins)[0], align='center', color='grey')
ax1.set_ylabel("% of realizations", fontsize=16)
ax1.set_xlabel("Years in decadal drought", fontsize=16)

ax2.barh(['Multi-decadal', 'Decadal'], [probability_multi, probability_decadal], height=0.2, align='center', color='grey')
ax2.set_xlabel("Fraction of all 100-year periods (%)", fontsize=16)

plt.show()

'''Extreme Value Analysis'''

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot()
ax.plot(Paleo['Year'][:429], Paleo['ScaledReconCisco'][:429], linewidth=3, color='black', label='Reconstructed paleo')
ax.plot(np.arange(1909, 2014), annual_historic_flows, linewidth=3, color='orange', label='Observed')

ax.hlines(y=np.min(annual_historic_flows),
          xmin=np.min(Paleo['Year'][:429]), xmax=2013,
          linewidth=2, linestyle='-.', color='orange')
ax.hlines(y=np.max(annual_historic_flows),
          xmin=np.min(Paleo['Year'][:429]), xmax=2013,
          linewidth=2, linestyle='-.', color='orange')
ax.hlines(y=mean_flow, xmin=np.min(Paleo['Year'][:429]), xmax=2013, linewidth=3, color='maroon', label='Mean')
ax.hlines(y=drought_threshold, xmin=np.min(Paleo['Year'][:429]), xmax=2013, linewidth=3,
          linestyle='--', color='maroon', label='Drought threshold')
ax.set_ylabel("Annual flow (Million $m^3$)", fontsize=16)
ax.set_xlabel("Year", fontsize=16)
ax.legend(fontsize=16, loc='upper left')
ax.tick_params(axis='both', labelsize=14)
plt.show()

# Frequency of deficits
drought_deficit = (drought_threshold-Paleo['ScaledReconCisco'][:429])/drought_threshold

hist, bins = np.histogram(drought_deficit, bins=30, density=True)
plt.bar(bins[:-1], hist, width=np.diff(bins)[0], align='edge', color='grey')
# mark 2002 deficit
deficit_2002 = (drought_threshold-annual_historic_flows[93])/drought_threshold
plt.axvline(x=deficit_2002, ymin=0, ymax=1, linewidth=3, linestyle='--', color='red', label='2002 deficit')
plt.ylabel("Frequency", fontsize=16)
plt.xlabel("Deficit ratio", fontsize=16)
plt.show()

# Distribution of extrema
#inverse_values = -1*Paleo['ScaledReconCisco'][:429]

block = 5
maxima = [max(drought_deficit[i:i+block]) for i in range(0, len(drought_deficit), block)]

fit = gev.fit(maxima)
shape, loc, scale = fit
min_interval, max_interval = gev.interval(0.99, shape, loc, scale)
xx = np.linspace(min_interval, max_interval, num=200)

yy = gev.pdf(xx, *fit)

hist, bins = np.histogram(maxima, bins=15, density=True)
plt.bar(bins[:-1], hist, width=np.diff(bins)[0], align='edge', color='grey')
plt.plot(xx, yy, linewidth=5, color='dodgerblue')
plt.axvline(x=deficit_2002, ymin=0, ymax=1, linewidth=3, linestyle='--', color='red', label='2002 deficit')
plt.ylabel("Probability density", fontsize=16)
plt.xlabel("Drought deficit ratio", fontsize=16)
plt.show()

# Calculate return periods for all extrema
maxima_sorted = np.sort(maxima)