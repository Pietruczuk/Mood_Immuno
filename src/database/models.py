from sqlalchemy import create_engine, Column, Integer, String, Text, Table, ForeignKey, Float, Date     # Importing basic SQLAlchemy components
from sqlalchemy.orm import declarative_base
from sqlalchemy.orm import relationship


# ---- Creating the Base class, which is required by SQLAlchemy to track all table definitions ------------------------/
Base = declarative_base()


# ---- Many-to-many relational table between articles and diseases ----------------------------------------------------/
article_disease_link = Table(
    'article_disease_link', Base.metadata,
    Column('article_id', Integer, ForeignKey('article.id'), primary_key=True),
    Column('disease_id', Integer, ForeignKey('disease.id'), primary_key=True)
)


# ---- Class representing the 'article' table in the database ---------------------------------------------------------/
class Article(Base):
    __tablename__ = 'article'

    id = Column(Integer, primary_key=True)     # Unique identifier in the table (INTEGER PRIMARY KEY)
    pubmed_id = Column(String, unique=True)    # PubMed ID of the article – unique for each entry in PubMed
    title = Column(Text)                       # Bibliographic data: title, abstract, journal, year of publication
    abstract = Column(Text)
    year = Column(Integer)
    doi = Column(String, unique=True)
    journal = Column(String)
    authors = Column(Text)                     # List of authors as text
    mesh_terms = Column(Text)                  # List of MeSH terms as text

    # Many-to-many relationship with 'disease', use of an intermediate table, synchronization
    diseases = relationship("Disease", secondary=article_disease_link, back_populates="articles")
    cytokines_entries = relationship("Cytokines", back_populates="article")
    cells_phenotypes_entries = relationship("CellsPhenotypes", back_populates="article")
    biochemistry_entries = relationship("Biochemistry", back_populates="article")
    patient_condition_entries = relationship("PatientCondition", back_populates="article")


# ---- Class representing the 'disease' table in the database ---------------------------------------------------------/
class Disease(Base):
    __tablename__ = 'disease'

    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True)

    # Many-to-many relationship with Article via an intermediate table
    articles = relationship("Article", secondary=article_disease_link, back_populates="diseases")

    # Relationships with other tables
    cells_phenotypes_entries = relationship("CellsPhenotypes", back_populates="disease")
    biochemistry_entries = relationship("Biochemistry", back_populates="disease")
    cytokines_entries = relationship("Cytokines", back_populates="disease")
    patient_condition_entries = relationship("PatientCondition", back_populates="disease")


# ---- Class representing the 'cytokines' table in the database -------------------------------------------------------/
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

    # Links to article
    article_id = Column(Integer, ForeignKey('article.id'))
    article = relationship("Article", back_populates="cytokines_entries")

    # Links to disease
    disease_id = Column(Integer, ForeignKey('disease.id'))
    disease = relationship("Disease", back_populates="cytokines_entries")

    # Links to cells_phenotypes
    cells_phenotypes_id = Column(Integer, ForeignKey('cells_phenotypes.id'))
    cells_phenotypes = relationship("CellsPhenotypes", back_populates="cytokines_entries")

    # Links to biochemistry
    biochemistry_id = Column(Integer, ForeignKey('biochemistry.id'))
    biochemistry = relationship("Biochemistry", back_populates="cytokines_entries")

    patient_condition_entries = relationship("PatientCondition", back_populates="cytokines")


# ---- Class representing the 'cells_phenotypes' table in the database ------------------------------------------------/
class CellsPhenotypes(Base):
    __tablename__ = 'cells_phenotypes'

    id = Column(Integer, primary_key=True)

    CD3CD4 = Column(String)
    CD3CD8 = Column(String)

    # Links to article
    article_id = Column(Integer, ForeignKey('article.id'))
    article = relationship("Article", back_populates="cells_phenotypes_entries")

    # Links to disease
    disease_id = Column(Integer, ForeignKey('disease.id'))
    disease = relationship("Disease", back_populates="cells_phenotypes_entries")

    # Links to cytokines
    cytokines_id = Column(Integer, ForeignKey('cytokines.id'))
    cytokines = relationship("Cytokines", back_populates="cells_phenotypes_entries")

    # Links to biochemistry
    biochemistry_id = Column(Integer, ForeignKey('biochemistry.id'))
    biochemistry = relationship("Biochemistry", back_populates="cells_phenotypes_entries")

    patient_condition_entries = relationship("PatientCondition", back_populates="cells_phenotypes")


# ---- Class representing the 'biochemistry' table in the database ----------------------------------------------------/
class Biochemistry(Base):
    __tablename__ = 'biochemistry'

    id = Column(Integer, primary_key=True)

    CRP = Column(String)

    # Links to article
    article_id = Column(Integer, ForeignKey('article.id'))
    article = relationship("Article", back_populates="biochemistry_entries")

    # Links to disease
    disease_id = Column(Integer, ForeignKey('disease.id'))
    disease = relationship("Disease", back_populates="biochemistry_entries")

    # Links to cytokines
    cytokines_id = Column(Integer, ForeignKey('cytokines.id'))
    cytokines = relationship("Cytokines", back_populates="biochemistry_entries")

    # Links to cells_phenotypes
    cells_phenotypes_id = Column(Integer, ForeignKey('cells_phenotypes.id'))
    cells_phenotypes = relationship("CellsPhenotypes", back_populates="biochemistry_entries")

    patient_condition_entries = relationship("PatientCondition", back_populates="biochemistry")



# ---- Class representing the 'patient_condition' table in the database -----------------------------------------------/
class PatientCondition(Base):
    __tablename__ = 'patient_condition'

    id = Column(Integer, primary_key=True)

    BDI = Column(Integer)      # Beck Depression Inventory
    HAM_D = Column(Integer)    # Hamilton Depression Rating Scale
    MADRS = Column(Integer)    # Montgomery–Åsberg Depression Rating Scale

    # Relationships with other tables
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
