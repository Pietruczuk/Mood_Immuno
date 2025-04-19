from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from src.database.models import Base

# ---- Funkcja tworząca połączenie z bazą ------------------------------------------------------------------------------
def init_db(db_path='c:/Users/amund/Desktop/github/Mood_Immuno/data/pubmed.db'):
    engine = create_engine(db_path)        # Tworzy silnik SQLAlchemy, który reprezentuje fizyczne połączenie z bazą danych
    Base.metadata.create_all(engine)       # Tworzy wszystkie tabele w bazie na podstawie klas zdefiniowanych w models.py, jeśli jeszcze nie istnieją.
    Session = sessionmaker(bind=engine)    # Tworzy klasę Session, która będzie używać naszego silnika jako źródła danych
    return Session()                       # Zwraca nową sesję (czyli połączenie z bazą) gotowe do użycia
