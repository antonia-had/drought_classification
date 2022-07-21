import numpy as np
import math
from matplotlib import pyplot as plt

def drought_counts(timeseries, drought_threshold, decadal_window=11, multidecadal_window=35):
    drought_windows = [multidecadal_window, decadal_window]
    all_years = len(timeseries)

    # Create list to store drought totals
    all_drought_years = []
    drought_instances = []

    for w in range(len(drought_windows)):
        window_size = drought_windows[w]
        # Calculate rolling means for each drought window
        i = 0
        # Initialize an empty list to store moving averages
        moving_averages = []

        # Loop through the array to consider
        # every window
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
    # length of multidecadal set and then length of decadal set
    total_drought_years = [len(years_list) for years_list in all_drought_years]
    return (total_drought_years)

def drought_identification_plots(timeseries, mean_flow, drought_threshold, decadal_window=11, multidecadal_window=35):
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
        # every window
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
    ax.set_xlabel("Year", fontsize=16)
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