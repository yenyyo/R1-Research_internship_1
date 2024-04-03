import matplotlib as plt

def pieplot(field):
    # Plot the pie chart
    plt.figure(figsize=(6, 6))
    plt.pie(field, labels=field.index, autopct='%1.1f%%', colors=['skyblue', 'lightcoral'])
    plt.title('Distribution of Sex')
    return plt

def barplot(series):
    # Create a bar chart
    plt.bar(series.index, series.values)

    # Set axis labels and title
    plt.xlabel("Racial Group")
    plt.ylabel("Count")
    plt.title("Racial Group Distribution")

    # Rotate x-axis labels for better readability with long labels
    plt.xticks(rotation=45)

    # Display the plot
    plt.show()