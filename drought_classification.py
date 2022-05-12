import numpy as np
import math
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import genextreme as gev

def drought_counts(timeseries, drought_threshold, decadal_window=11, multidecadal_window=35):
    drought_windows = [multidecadal_window, decadal_window]
    all_years = len(timeseries)

    # Create list to store drought totals
    all_drought_years = []

    for w in range(len(drought_windows)):
        window_size = drought_windows[w]
        # Calculate rolling means for each drought window
        i = 0
        # Initialize an empty list to store moving averages
        moving_averages = []

        # Loop through the array to consider
        # every window of size 3
        while i < len(timeseries) - window_size + 1:
            # Store elements from i to i+window_size
            # in list to get the current window
            window = timeseries[i: i + window_size]
            # Calculate the average of current window
            window_average = round(sum(window) / window_size, 2)
            # Store the average of current
            # window in moving average list
            moving_averages.append(window_average)
            # Shift window to right by one position
            i += 1

        # Plot results
        half_window = math.floor(window_size / 2)

        highlight_years = np.arange(half_window, all_years - half_window)[moving_averages < drought_threshold]
        all_drought_years_in_window = []
        for year in highlight_years:
            first_year = year - half_window
            last_year = year + half_window
            all_drought_years_in_window.extend(np.arange(first_year, last_year + 1))

        all_drought_years.append(set(all_drought_years_in_window))
    # lenght of multidecadal set and then length of decadal set
    total_drought_years = [len(years_list) for years_list in all_drought_years]
    return (total_drought_years)

def drought_identification_plots(timeseries, drought_threshold, decadal_window=11, multidecadal_window=35):
    drought_windows = [multidecadal_window, decadal_window]
    drought_colors = ['#CA6702', '#EE9B00']
    all_years = len(timeseries)

    fig = plt.figure(figsize=(16, 9))
    # parameters to specify the width and height ratios between rows and columns
    widths = [1]
    heights = [5, 1]

    gspec = fig.add_gridspec(ncols=1, nrows=2, width_ratios=widths, height_ratios=heights)
    ax = fig.add_subplot(gspec[0, 0])
    ax.plot(np.arange(all_years), timeseries, linewidth=3, color='#001219')
    ax.hlines(y=mean_flow, xmin=0, xmax=len(timeseries), linewidth=2, color='#005F73', label='Mean')
    ax.hlines(y=drought_threshold, xmin=0, xmax=len(timeseries), linewidth=2,
              linestyle='--', color='#005F73', label='Drought threshold')

    # Create list to store drought totals
    all_drought_years = []

    for w in range(len(drought_windows)):
        window_size = drought_windows[w]
        # Calculate rolling means for each drought window
        i = 0
        # Initialize an empty list to store moving averages
        moving_averages = []

        # Loop through the array to consider
        # every window of size 3
        while i < len(timeseries) - window_size + 1:
            # Store elements from i to i+window_size
            # in list to get the current window
            window = timeseries[i: i + window_size]

            # Calculate the average of current window
            window_average = round(sum(window) / window_size, 2)

            # Store the average of current
            # window in moving average list
            moving_averages.append(window_average)

            # Shift window to right by one position
            i += 1

        # Plot results
        half_window = math.floor(window_size / 2)

        highlight_years = np.arange(half_window, all_years - half_window)[moving_averages < drought_threshold]
        all_drought_years_in_window = []
        for year in highlight_years:
            first_year = year - half_window
            last_year = year + half_window
            # the window to highlight is +/- 5 years from the identified year crossing the threshold
            ax.axvspan(first_year, last_year, color=drought_colors[w])
            #         ax.plot(np.arange(half_window, all_years-half_window),
            #                 moving_averages, linewidth=3, color='grey', label=f'{window_size}-yr mean')
            all_drought_years_in_window.extend(np.arange(first_year, last_year + 1))

        all_drought_years.append(set(all_drought_years_in_window))
    # remove decadal drought years that also appear in multidecadal
    all_drought_years[1] = set(all_drought_years[1]) - set(all_drought_years[0])
    # lenght of multidecadal set and then length of decadal set
    total_drought_years = [len(years_list) for years_list in all_drought_years]
    ax.tick_params(axis='both', labelsize=14)
    ax.set_ylabel("Annual flow (Million $m^3$)", fontsize=16)
    ax.set_xlabel("Year in realization", fontsize=16)
    ax.legend(fontsize=16, loc='upper right')

    ax2 = fig.add_subplot(gspec[1, :])
    # plot decadal first and then multidecadal
    ax2.barh(1, total_drought_years[1], color='#EE9B00')
    ax2.barh(1, total_drought_years[0], left=total_drought_years[1], color='#CA6702')

    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_yticks([])
    ax2.tick_params(axis='x', labelsize=14)
    ax2.set_xlabel("Total years in drought", fontsize=16)
    ax2.set_xlim(ax.get_xlim())
    plt.show()
    return (total_drought_years)


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

drought_identification_plots(Paleo['ScaledReconCisco'][:429], drought_threshold)

# Repeat for synthetically generated flows
paleo_flows = np.load('./data/Paleo_SOWs_flows.npy')
paleo_flows_annual = np.sum(paleo_flows, axis=2)*1233.4818/1000000

multi_numbers = []
decadal_numbers = []
for j in range(len(paleo_flows_annual[:,0])):
    [multi, decadal] = drought_counts(paleo_flows_annual[j,:], drought_threshold, decadal_window=11, multidecadal_window=35)
    multi_numbers.append(multi)
    decadal_numbers.append(decadal)

zero_droughts=[i for i, e in enumerate(decadal_numbers) if e == 0]
probability_decadal = 100*(1-len(zero_droughts)/len(decadal_numbers))

zero_multi=[i for i, e in enumerate(multi_numbers) if e == 0]
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