import requests

url = "https://dfam.org/api/families"
params = {
    # The summary format is metadata-only and does not include
    # full details such as the consensus sequence and citations
    "format": "summary",

    # Only retrieve the first 10 results in this query
    "limit": "10",

    # Search in Caenorhabditis elegans (worm)
    "clade": 96803,

    # Include families from ancestor and descendant taxa in the results
    "clade_relatives": "both",
}

response = requests.get(url, params=params)
results = response.json()["results"]

# Prints "Vingi-2_CE" at the time of this writing
print(results)

url = "https://dfam.org/api/classes"
response = requests.get(url)
results = response.json()
print(results)
