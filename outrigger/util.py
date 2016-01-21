import datetime

def timestamp():
    return str(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
