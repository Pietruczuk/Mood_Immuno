from sqlalchemy import create_engine, Column, Integer, String, Text, Table, ForeignKey, Float, Date     # Import podstawowych komponentów SQLAlchemy
from sqlalchemy.orm import declarative_base
from sqlalchemy.orm import relationship


# ---- Tworzenie klasy Base, która jest wymagana przez SQLAlchemy do śledzenia wszystkich definicji tabel -------------
Base = declarative_base()


# ---- Tabela relacyjna wiele-do-wielu między artykułami a chorobami --------------------------------------------------
article_disease_link = Table(
    'article_disease_link', Base.metadata,
    Column('article_id', Integer, ForeignKey('article.id'), primary_key=True),
    Column('disease_id', Integer, ForeignKey('disease.id'), primary_key=True)
)


# ---- Klasa reprezentująca tabelę 'article' w bazie danych -----------------------------------------------------------
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
    immuno_entries = relationship("Immuno", back_populates="article")


# ---- Klasa reprezentująca tabelę 'disease' w bazie danych -----------------------------------------------------------
class Disease(Base):
    __tablename__ = 'disease'

    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True)

    articles = relationship("Article", secondary=article_disease_link, back_populates="diseases")


# ---- Klasa reprezentująca tabelę 'immuno' w bazie danych -----------------------------------------------------------
class Immuno(Base):
    __tablename__ = 'immuno'

    id = Column(Integer, primary_key=True)
    nazwa = Column(String)
    CRP = Column(Float)
    IL1 = Column(Float)
    IL6 = Column(Float)

    article_id = Column(Integer, ForeignKey('article.id'))  # klucz obcy do article
    article = relationship("Article", back_populates="immuno_entries")  # relacja odwrotna