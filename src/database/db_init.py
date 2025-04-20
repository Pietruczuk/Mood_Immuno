from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from src.database.models import Base

# ---- Function creating a connection to the database -----------------------------------------------------------------/
def init_db(db_path='c:/Users/amund/Desktop/github/Mood_Immuno/data/pubmed.db'):
    engine = create_engine(db_path)        # Creates a SQLAlchemy engine that represents a physical connection to the database.
    Base.metadata.create_all(engine)       # Creates all tables in the database based on classes defined in models.py, if they do not already exist.
    Session = sessionmaker(bind=engine)    # Creates a Session class that will use our engine as a data source
    return Session()                       # Returns a new session (i.e., a connection to the database) ready for use.
