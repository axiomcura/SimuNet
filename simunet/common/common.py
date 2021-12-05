from datetime import datetime

def generate_id() -> str :
    """Generates a unique id based on time and date
    """
    unique_id = datetime.now().strftime("%M%d%y-%H%M%S")
    return unique_id
    