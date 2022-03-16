
from binance.client import Client as BiCl

class CoinTrader:

    def __init__(self):
        pass


if __name__ == '__main__':
    # TODO
    # --- connect to binance api ---
    API_KEY = "TKZaR8z9wsbfQ6MQHTh7HCDJ70mezWcd8Ie0RuaYOPY2cwB7zMIKiDxrQw2pTAGJ"
    API_SECRET = "QBosTWqegkBdWJJ1XDdD4k3Qx7Bn8MxSFQzvQXZybHxJraUqiyF8AG0utrD7J2Wb"

    client = BiCl(api_key=API_KEY,api_secret=API_SECRET)
    stat = client.get_account_status()
    acc = client.get_account()
    print(acc) # OK
    print(stat) # OK

    # get coinprices
    info = client.get_symbol_info("BTCUSDT")
    print(info)
    #ticker = client.get_symbol_ticker()
    #print(ticker)
    ex_info = client.get_exchange_info()
    #print(ex_info)

    # get coinprices for fluctuations

    # get recommendations based on 1 day , 1 week , 1 month fluctuations !

    # apply trading strategy ..
    # trade...
    # save trading logs

    # compare results to other strategies

