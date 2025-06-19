###We plot scaling analysis results for HIV Drug Resistance tool Snakemake pipeline
import matplotlib.pyplot as plt
import pandas as pd

scaling_data = pd.read_csv('scaling_analysis_data.tsv', sep='\t')

# Convert time strings to seconds for plotting
def convert_time_to_seconds(time_str):
    if 'h' in time_str:
        hours, minutes = time_str.split('h')
        minutes, seconds = minutes.split('m')
        seconds = seconds.replace('s', '').strip()
        total_seconds = int(hours) * 3600 + int(minutes) * 60 + int(seconds)
    else:
        if 'm' in time_str:
            minutes, seconds = time_str.split('m')
            seconds = seconds.replace('s', '').strip()
            total_seconds = int(minutes) * 60 + int(seconds)
        else:
            total_seconds = int(time_str.replace('s', '').strip())
    return total_seconds

# Convert seconds back to a readable format for y-axis labels
def convert_seconds_to_time_format(seconds):
    if seconds >= 3600:
        hours = seconds // 3600
        seconds %= 3600
        minutes = seconds // 60
        seconds %= 60
        return f"{hours}h {minutes}m"
    elif seconds >= 60:
        minutes = seconds // 60
        seconds %= 60
        return f"{minutes}m"
    else:
        return f"{seconds}s"

scaling_data['time_4_cores'] = scaling_data['time_4_cores'].apply(convert_time_to_seconds)
scaling_data['time_32_cores'] = scaling_data['time_32_cores'].apply(convert_time_to_seconds)

##We plot the scaling analysis results (Figure 4)
plt.figure(figsize=(10, 6))
plt.plot(scaling_data['n_samples'], scaling_data['time_4_cores'], marker='o', label='4 cores')
plt.plot(scaling_data['n_samples'], scaling_data['time_32_cores'], marker='o', label='32 cores')
plt.xlabel('Number of samples')
plt.ylabel('Time')
plt.title('Scaling Analysis of HIV Drug Resistance Tool')
plt.xticks([1, 250, 500, 750, 1000])
yticks = [0, 1800, 3600, 5400]  # in seconds
plt.yticks(yticks, [convert_seconds_to_time_format(tick) for tick in yticks])
plt.legend()
plt.savefig('figures/scaling_analysis_plot.png', dpi=300)
