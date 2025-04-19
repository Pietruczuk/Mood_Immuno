from sqlalchemy import create_engine, Column, Integer, String, Text, Table, ForeignKey, Float, Date     # Import podstawowych komponentów SQLAlchemy
from sqlalchemy.orm import declarative_base
from sqlalchemy.orm import relationship


# ---- Tworzenie klasy Base, która jest wymagana przez SQLAlchemy do śledzenia wszystkich definicji tabel -------------/
Base = declarative_base()


# ---- Tabela relacyjna wiele-do-wielu między artykułami a chorobami --------------------------------------------------/
article_disease_link = Table(
    'article_disease_link', Base.metadata,
    Column('article_id', Integer, ForeignKey('article.id'), primary_key=True),
    Column('disease_id', Integer, ForeignKey('disease.id'), primary_key=True)
)


# ---- Klasa reprezentująca tabelę 'article' w bazie danych -----------------------------------------------------------/
class Article(Base):
    __tablename__ = 'article'

    id = Column(Integer, primary_key=True)     # Unikalny identyfikator w tabeli (INTEGER PRIMARY KEY)
    pubmed_id = Column(String, unique=True)    # PubMed ID artykułu – unikalny dla każdego wpisu w Pubmed
    title = Column(Text)                       # Dane bibliograficzne: tytuł, streszczenie, czasopismo, rok publikacji
    abstract = Column(Text)
    year = Column(Integer)
    doi = Column(String, unique=True)
    journal = Column(String)
    authors = Column(Text)                     # Lista autorów jako tekst
    mesh_terms = Column(Text)                  # Lista MeSH terms jako tekst

    # Relacja wiele-do-wielu z 'disease', użycie tabeli pośredniczącej, synchronizacja
    diseases = relationship("Disease", secondary=article_disease_link, back_populates="articles")
    cytokines_entries = relationship("Cytokines", back_populates="article")
    cells_phenotypes_entries = relationship("CellsPhenotypes", back_populates="article")
    biochemistry_entries = relationship("Biochemistry", back_populates="article")
    patient_condition_entries = relationship("PatientCondition", back_populates="article")


# ---- Klasa reprezentująca tabelę 'disease' w bazie danych -----------------------------------------------------------/
class Disease(Base):
    __tablename__ = 'disease'

    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True)

    # Relacja wiele-do-wielu z Article przez tabelę pośrednią
    articles = relationship("Article", secondary=article_disease_link, back_populates="diseases")

    # Relacje z innymi tabelami
    cells_phenotypes_entries = relationship("CellsPhenotypes", back_populates="disease")
    biochemistry_entries = relationship("Biochemistry", back_populates="disease")
    cytokines_entries = relationship("Cytokines", back_populates="disease")
    patient_condition_entries = relationship("PatientCondition", back_populates="disease")


# ---- Klasa reprezentująca tabelę 'cytokines' w bazie danych ---------------------------------------------------------/
class Cytokines(Base):
    __tablename__ = 'cytokines'

    id = Column(Integer, primary_key=True)

    IL1 = Column(String)
    IL6 = Column(String)
    IL8 = Column(String)
    IL9 = Column(String)
    IL12 = Column(String)
    IL15 = Column(String)
    IL17 = Column(String)
    IL18 = Column(String)
    TNF_alpha = Column(String)
    TNF_beta = Column(String)
    INF_gamma = Column(String)

    # Powiązania z article
    article_id = Column(Integer, ForeignKey('article.id'))
    article = relationship("Article", back_populates="cytokines_entries")

    # Powiązania z disease
    disease_id = Column(Integer, ForeignKey('disease.id'))
    disease = relationship("Disease", back_populates="cytokines_entries")

    # Powiązania z cells_phenotypes
    cells_phenotypes_id = Column(Integer, ForeignKey('cells_phenotypes.id'))
    cells_phenotypes = relationship("CellsPhenotypes", back_populates="cytokines_entries")

    # Powiązania z biochemistry
    biochemistry_id = Column(Integer, ForeignKey('biochemistry.id'))
    biochemistry = relationship("Biochemistry", back_populates="cytokines_entries")

    patient_condition_entries = relationship("PatientCondition", back_populates="cytokines")


# ---- Klasa reprezentująca tabelę 'cells_phenotypes' w bazie danych --------------------------------------------------/
class CellsPhenotypes(Base):
    __tablename__ = 'cells_phenotypes'

    id = Column(Integer, primary_key=True)

    CD3CD4 = Column(String)
    CD3CD8 = Column(String)

    # Powiązania z tabelą article
    article_id = Column(Integer, ForeignKey('article.id'))
    article = relationship("Article", back_populates="cells_phenotypes_entries")

    # Powiązania z tabelą disease
    disease_id = Column(Integer, ForeignKey('disease.id'))
    disease = relationship("Disease", back_populates="cells_phenotypes_entries")

    # Powiązania z tabelą cytokines
    cytokines_id = Column(Integer, ForeignKey('cytokines.id'))
    cytokines = relationship("Cytokines", back_populates="cells_phenotypes_entries")

    # Powiązania z tabelą biochemistry
    biochemistry_id = Column(Integer, ForeignKey('biochemistry.id'))
    biochemistry = relationship("Biochemistry", back_populates="cells_phenotypes_entries")

    patient_condition_entries = relationship("PatientCondition", back_populates="cells_phenotypes")


# ---- Klasa reprezentująca tabelę 'biochemistry' w bazie danych --------------------------------------------------/
class Biochemistry(Base):
    __tablename__ = 'biochemistry'

    id = Column(Integer, primary_key=True)

    CRP = Column(String)

    # Powiązania z article
    article_id = Column(Integer, ForeignKey('article.id'))
    article = relationship("Article", back_populates="biochemistry_entries")

    # Powiązania z disease
    disease_id = Column(Integer, ForeignKey('disease.id'))
    disease = relationship("Disease", back_populates="biochemistry_entries")

    # Powiązania z cytokines
    cytokines_id = Column(Integer, ForeignKey('cytokines.id'))
    cytokines = relationship("Cytokines", back_populates="biochemistry_entries")

    # Powiązania z cells_phenotypes
    cells_phenotypes_id = Column(Integer, ForeignKey('cells_phenotypes.id'))
    cells_phenotypes = relationship("CellsPhenotypes", back_populates="biochemistry_entries")

    patient_condition_entries = relationship("PatientCondition", back_populates="biochemistry")



# ---- Klasa reprezentująca tabelę 'patient_condition' w bazie danych -------------------------------------------------/
class PatientCondition(Base):
    __tablename__ = 'patient_condition'

    id = Column(Integer, primary_key=True)

    BDI = Column(Integer)      # Beck Depression Inventory
    HAM_D = Column(Integer)    # Hamilton Depression Rating Scale
    MADRS = Column(Integer)    # Montgomery–Åsberg Depression Rating Scale

    # Relacje z pozostałymi tabelami:
    article_id = Column(Integer, ForeignKey('article.id'))
    article = relationship("Article", back_populates="patient_condition_entries")

    disease_id = Column(Integer, ForeignKey('disease.id'))
    disease = relationship("Disease", back_populates="patient_condition_entries")

    cytokines_id = Column(Integer, ForeignKey('cytokines.id'))
    cytokines = relationship("Cytokines", back_populates="patient_condition_entries")

    cells_phenotypes_id = Column(Integer, ForeignKey('cells_phenotypes.id'))
    cells_phenotypes = relationship("CellsPhenotypes", back_populates="patient_condition_entries")

    biochemistry_id = Column(Integer, ForeignKey('biochemistry.id'))
    biochemistry = relationship("Biochemistry", back_populates="patient_condition_entries")
