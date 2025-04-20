# ---- Code for downloading data from Pubmed --------------------------------------------------------------------------/
# ---- Autor: Krzysztof Pietruczuk ------------------------------------------------------------------------------------/

from Bio import Entrez                   # A module from the Biopython library enabling communication with NCBI servers.


# ---- Configuration: provide your email address so that NCBI knows who is using the API ------------------------------/
Entrez.email = "krzysztof.pietruczuk@gumed.edu.pl"


# ---- Function for searching articles in PubMed based on a query -----------------------------------------------------/
def fetch_pubmed_ids(query, max_results=100):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)          # Submitting a query to PubMed
    record = Entrez.read(handle)               # Parses the result in XML format and converts it to a Python data structure.
    return record["IdList"]                    # Returns a list of PubMed IDs – unique identifiers for articles


# ---- Function that retrieves detailed data based on PubMed IDs of articles ------------------------------------------/
def fetch_details(pubmed_id: str):
    handle = Entrez.efetch(db="pubmed", id=pubmed_id, rettype="xml")
    record = Entrez.read(handle)

    # Check if this is a PubmedArticle from MedlineCitation
    if isinstance(record, dict) and "PubmedArticle" in record:
        record = record["PubmedArticle"][0]
    else:
        raise ValueError(f"⚠ No data from NCBI for PubMed ID: {pubmed_id}. Answer: {record}")

    citation = record.get("MedlineCitation", {})
    article_info = citation.get("Article", {})

    title = article_info.get("ArticleTitle", "No title")
    abstract_list = article_info.get("Abstract", {}).get("AbstractText", ["No abstract"])
    abstract = abstract_list[0] if isinstance(abstract_list, list) else abstract_list

    pub_date = article_info.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
    year_str = pub_date.get("Year", "")
    year = int(year_str) if year_str.isdigit() else None

    journal = article_info.get("Journal", {}).get("Title", "unknown")

    doi = None
    id_list = record.get("PubmedData", {}).get("ArticleIdList", [])
    for item in id_list:
        if hasattr(item, "attributes") and item.attributes.get("IdType") == "doi":
            doi = str(item)
            break

    authors = []
    for author in article_info.get("AuthorList", []):
        try:
            firstname = author.get("ForeName", "")
            lastname = author.get("LastName", "")
            full_name = f"{firstname} {lastname}".strip()
            if full_name:
                authors.append(full_name)
        except Exception:
            continue

    mesh_terms = []
    for mesh in citation.get("MeshHeadingList", []):
        descriptor = mesh.get("DescriptorName")
        if descriptor:
            mesh_terms.append(str(descriptor))

    return {
        "pubmed_id": pubmed_id,
        "title": title,
        "abstract": abstract,
        "year": year,
        "doi": doi,
        "journal": journal,
        "authors": ", ".join(authors),
        "mesh_terms": ", ".join(mesh_terms)
    }
