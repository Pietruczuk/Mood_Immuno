from sqlalchemy.exc import IntegrityError
from src.database.db_init import init_db
from src.database.models import Article, Disease
from src.pubmed.fetch_pubmed import fetch_pubmed_ids, fetch_details

import os
import traceback

# --- Konfiguracja zapytania ---
TERM = "depressive disorder"
QUERY = f'"{TERM}"[MeSH Terms] AND "humans"[MeSH Terms] AND ("2014"[dp] : "2024"[dp])'
MAX_RESULTS = 7

# --- Połączenie z bazą ---
db_path = os.path.abspath("data/pubmed.db")
session = init_db(f"sqlite:///{db_path}")
print(f"🗃 Połączono z bazą danych: {db_path}")

# --- Upewnij się, że choroba jest w bazie ---
disease = session.query(Disease).filter_by(name=TERM).first()
if not disease:
    disease = Disease(name=TERM)
    session.add(disease)
    session.commit()
    print(f"✅ Dodano nową chorobę: {TERM}")

# --- Pobieranie ID ---
ids = fetch_pubmed_ids(QUERY, max_results=MAX_RESULTS)
print(f"🔍 Pobrano {len(ids)} PubMed ID")

# --- Przetwarzanie artykułów ---
for pubmed_id in ids:
    if session.query(Article).filter_by(pubmed_id=pubmed_id).first():
        print(f"⚠ Artykuł {pubmed_id} już istnieje — pominięto.")
        continue

    try:
        data = fetch_details(pubmed_id)

        article = Article(
            pubmed_id=data["pubmed_id"],
            title=data["title"],
            abstract=data["abstract"],
            year=data["year"],
            doi=data["doi"],
            journal=data["journal"],
            authors=data["authors"],
            mesh_terms=data["mesh_terms"]
        )

        article.diseases.append(disease)
        session.add(article)
        session.commit()
        print(f"✅ Dodano artykuł: {pubmed_id}")

    except Exception as e:
        session.rollback()
        print(f"❌ Błąd przy dodawaniu artykułu {pubmed_id}:")
        traceback.print_exc()

print("🎉 Gotowe!")
