# ---- Kod pobierania danych z Pubmed -------
# ---- Autor: Krzysztof Pietruczuk ----------

from Bio import Entrez                        # Moduł z biblioteki Biopython umożliwiający komunikację z serwerami NCBI


# ---- Konfiguracja: podanie e-maila, aby NCBI wiedziało, kto korzysta z API ------------------------------------------
Entrez.email = "krzysztof.pietruczuk@gumed.edu.pl"


# ---- Funkcja wyszukująca artykuły w PubMed na podstawie zapytania ---------------------------------------------------
def fetch_pubmed_ids(query, max_results=100):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)  # Wysłanie zapytania do PubMed
    record = Entrez.read(handle)               # Parsuje wynik w formacie XML i zamienia go na strukturę danych Pythona
    return record["IdList"]                    # Zwraca listę PubMed ID – unikalnych identyfikatorów artykułów


# ---- Funkcja pobierająca szczegółowe dane na podstawie PubMed ID artykułów ------------------------------------------
def fetch_details(pubmed_id: str):
    handle = Entrez.efetch(db="pubmed", id=pubmed_id, rettype="xml")
    record = Entrez.read(handle)

    # Sprawdź, czy to PubmedArticle z MedlineCitation
    if isinstance(record, dict) and "PubmedArticle" in record:
        record = record["PubmedArticle"][0]
    else:
        raise ValueError(f"⚠ Brak danych z NCBI dla PubMed ID: {pubmed_id}. Odpowiedź: {record}")

    citation = record.get("MedlineCitation", {})
    article_info = citation.get("Article", {})

    title = article_info.get("ArticleTitle", "Brak tytułu")
    abstract_list = article_info.get("Abstract", {}).get("AbstractText", ["Brak abstraktu"])
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
