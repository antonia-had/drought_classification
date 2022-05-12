import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import stats
import pandas as pd
from matplotlib import pyplot as plt
import scipy.stats as ss


def ba_nonparametric_qm(sim, sim_ref, obs_ref, nbins=20):
    # from http://www.pik-potsdam.de/~menz/IKI-Oasis/capacity_building/scripts/
    order = 1
    bins = np.linspace(0, 100, nbins)

    bin_obs_ref = np.percentile(obs_ref, q=bins)
    bin_sim_ref = np.percentile(sim_ref, q=bins)

    g = InterpolatedUnivariateSpline(bin_sim_ref, bin_obs_ref, k=order)

    return g(sim)

def ba_parametric_qm(sim, sim_ref, obs_ref):
    # from http://www.pik-potsdam.de/~menz/IKI-Oasis/capacity_building/scripts/

    dist_obs_ref_mean, dist_obs_ref_std = stats.norm.fit(obs_ref)
    dist_sim_ref_mean, dist_sim_ref_std = stats.norm.fit(sim_ref)

    F_sim = stats.norm.cdf(sim, dist_sim_ref_mean, dist_sim_ref_std)

    return stats.norm.ppf(F_sim, dist_obs_ref_mean, dist_obs_ref_std)

# load paleo data at Cisco
Paleo = pd.read_csv('./data/Reconstruction/Cisco_Recon_v_Observed_v_Stateline.csv')

# re-scale Cisco data to estimate data at CO-UT state line
factor = np.nanmean(Paleo['ObservedNaturalStateline']/Paleo['ObservedNaturalCisco'])
Paleo['ScaledNaturalCisco'] = Paleo['ObservedNaturalCisco']*factor
Paleo['ScaledReconCisco'] = Paleo['ReconCisco']*factor

#compare reconstructed with observed for overlapping period
fig, (ax1) = plt.subplots(1,1, figsize=(12,8))
ax1.scatter(np.sort(Paleo['ScaledNaturalCisco'][337:426]), np.sort(Paleo['ScaledReconCisco'][337:426]), color='royalblue', alpha=1)
ax1.plot([np.min(Paleo['ScaledNaturalCisco'][337:426]), np.max(Paleo['ScaledNaturalCisco'][337:426])],
         [np.min(Paleo['ScaledNaturalCisco'][337:426]), np.max(Paleo['ScaledNaturalCisco'][337:426])], color='r')
ax1.set_xlabel('Observed flows at Cisco', fontsize=16)
ax1.set_ylabel('Reconstructed flows at Cisco', fontsize=16)
plt.show()

#compare reconstructed + noise with observed for overlapping period
nsims=100
stdev = np.std(Paleo['FractionScalingResid'][337:426])
fig, (ax1) = plt.subplots(1,1, figsize=(12,8))
for i in range(nsims):
    flows = Paleo['ScaledReconCisco'][337:426] + Paleo['ScaledReconCisco'][337:426]*ss.norm.rvs(0,stdev,426-337)
    ax1.scatter(np.sort(Paleo['ScaledNaturalCisco'][337:426]), np.sort(flows), color='royalblue', alpha=0.5)
ax1.plot([np.min(Paleo['ScaledNaturalCisco'][337:426]), np.max(Paleo['ScaledNaturalCisco'][337:426])],
         [np.min(Paleo['ScaledNaturalCisco'][337:426]), np.max(Paleo['ScaledNaturalCisco'][337:426])], color='r')
ax1.set_xlabel('Observed flows at Cisco', fontsize=16)
ax1.set_ylabel('Reconstructed flows at Cisco', fontsize=16)
plt.show()

#plot flow duration curves for paleo and observations (overlapping years)
P_paleo = np.arange(1.,len(Paleo['ScaledReconCisco'][337:426])+1)*100 / len(Paleo['ScaledReconCisco'][337:426])
P_observed = np.arange(1.,len(Paleo['ScaledNaturalCisco'][337:426])+1)*100 / len(Paleo['ScaledNaturalCisco'][337:426])

fig, (ax1) = plt.subplots(1,1, figsize=(12,8))
ax1.plot(P_paleo, np.sort(Paleo['ScaledReconCisco'][337:426])[::-1], linewidth=3, color='royalblue', label='Paleo at Cisco')
ax1.plot(P_observed, np.sort(Paleo['ScaledNaturalCisco'][337:426])[::-1], linewidth=3, color='red', label='Observed at Cisco')
ax1.set_ylabel('Flows (m$^3$)', fontsize=16)
ax1.set_xlabel('Excedence probability', fontsize=16)
plt.legend(fontsize=16)
plt.show()

# adjust paleo data with quantile mapping (using only overlap)
adjusted_overlap = ba_nonparametric_qm(Paleo['ScaledReconCisco'][337:426], Paleo['ScaledReconCisco'][337:426], Paleo['ScaledNaturalCisco'][337:426], nbins=20)

fig, (ax1) = plt.subplots(1,1, figsize=(12,8))
ax1.scatter(np.sort(Paleo['ScaledNaturalCisco'][337:426]), np.sort(adjusted_overlap), color='royalblue', alpha=1)
ax1.plot([np.min(Paleo['ScaledNaturalCisco'][337:426]), np.max(Paleo['ScaledNaturalCisco'][337:426])],
         [np.min(Paleo['ScaledNaturalCisco'][337:426]), np.max(Paleo['ScaledNaturalCisco'][337:426])], color='r')
ax1.set_xlabel('Observed flows at Cisco ($m^3$)', fontsize=16)
ax1.set_ylabel('Biased corrected flows at Cisco ($m^3$)', fontsize=16)
plt.show()

fig, (ax1) = plt.subplots(1,1, figsize=(12,8))
ax1.plot(P_paleo, np.sort(adjusted_overlap)[::-1], linewidth=3, color='royalblue', label='Adjusted paleo at Cisco')
ax1.plot(P_observed, np.sort(Paleo['ScaledNaturalCisco'][337:426])[::-1], linewidth=3, color='red', label='Observed at Cisco')
ax1.set_ylabel('Flows (m$^3$)', fontsize=16)
ax1.set_xlabel('Excedence probability', fontsize=16)
plt.legend(fontsize=16)
plt.show()

# adjust paleo data with quantile mapping (using all the record)

adjusted = ba_parametric_qm(Paleo['ScaledReconCisco'][:429], Paleo['ScaledReconCisco'][337:426], Paleo['ScaledNaturalCisco'][337:426])


P_adjusted = np.arange(1., len(adjusted)+1)*100 / len(adjusted)

fig, (ax1) = plt.subplots(1,1, figsize=(12,8))
ax1.plot(P_adjusted, np.sort(adjusted)[::-1], linewidth=3, color='orange', label='Adjusted paleo at Cisco (full record)')
ax1.plot(P_adjusted, np.sort(Paleo['ScaledReconCisco'][:429])[::-1], linewidth=3, color='green', label='Original paleo at Cisco (full record)')
#ax1.plot(P_observed, np.sort(Paleo['ScaledNaturalCisco'][337:426])[::-1], linewidth=3, color='red', label='Observed at Cisco')
ax1.set_ylabel('Flows (m$^3$)', fontsize=16)
ax1.set_xlabel('Excedence probability', fontsize=16)
plt.legend(fontsize=16)
plt.show()