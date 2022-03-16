# SOURCES :
#https://python-binance.readthedocs.io/en/latest/
#https://medium.com/swlh/retrieving-full-historical-data-for-every-cryptocurrency-on-binance-bitmex-using-the-python-apis-27b47fd8137f

from datetime import timedelta, datetime
import math
import os
import matplotlib.pyplot as plt
import parser
import numpy as np
import requests
import pandas as pd
import json
from binance.client import Client

def decideBuy(ref:float, actual:float, capital : float = 1000):

    # buy or not ?
    # 10% down criterium
    buy : bool = False
    threshold = 0.05 # 10% down
    if actual/ref < (1-threshold):
        buy = True
    return buy

def decideSell(ref:float, actual:float, btc : float = 1000):

    # sell or not ?
    # 25% up criterium
    sell : bool = False
    threshold = 0.1 # 10% up
    if actual/ref > (1+threshold):
        sell = True
    return sell

#copied from SOURCES
def get_all_binance(symbol, kline_size, save = False):
    filename = '%s-%s-data.csv' % (symbol, kline_size)
    if os.path.isfile(filename): data_df = pd.read_csv(filename)
    else: data_df = pd.DataFrame()
    oldest_point, newest_point = minutes_of_new_data(symbol, kline_size, data_df, source = "binance")
    delta_min = (newest_point - oldest_point).total_seconds()/60
    available_data = math.ceil(delta_min/_binsizes[kline_size])
    if oldest_point == datetime.strptime('1 Jan 2017', '%d %b %Y'): print('Downloading all available %s data for %s. Be patient..!' % (kline_size, symbol))
    else: print('Downloading %d minutes of new data available for %s, i.e. %d instances of %s data.' % (delta_min, symbol, available_data, kline_size))
    klines = binance_client.get_historical_klines(symbol, kline_size, oldest_point.strftime("%d %b %Y %H:%M:%S"), newest_point.strftime("%d %b %Y %H:%M:%S"))
    data = pd.DataFrame(klines, columns = ['timestamp', 'open', 'high', 'low', 'close', 'volume', 'close_time', 'quote_av', 'trades', 'tb_base_av', 'tb_quote_av', 'ignore' ])
    data['timestamp'] = pd.to_datetime(data['timestamp'], unit='ms')
    if len(data_df) > 0:
        temp_df = pd.DataFrame(data)
        data_df = data_df.append(temp_df)
    else: data_df = data
    data_df.set_index('timestamp', inplace=True)
    if save: data_df.to_csv(filename)
    print('All caught up..!')
    return data_df

#copied from SOURCES
def minutes_of_new_data(symbol, kline_size, data, source):
    if len(data) > 0:  old = parser.parse(data["timestamp"].iloc[-1])
    elif source == "binance": old = datetime.strptime('1 Jan 2017', '%d %b %Y')
    if source == "binance": new = pd.to_datetime(binance_client.get_klines(symbol=symbol, interval=kline_size)[-1][0], unit='ms')
    return old, new

def derivative_array(arr, delta = 2):
    der = arr.copy()
    for i in range(0,len(der)):
        if i-delta <= 0 or i+delta >= len(der):
            der[i] = -1
        else:
            der[i] = der[i+delta] - der[i-delta]
    return der

# get momentarily price of CRYPTO
def get_latest_price(crypto):
    session = requests.Session()
    session.trust_env = False;
    response = session.get(_TICKER_API_URL+crypto)
    response_json = response.json()
    return response_json['price']

# get average price of last n_days of price_arr
def get_averages(n_days,price_arr):
    avg_ndays = price_arr.copy();
    for i in range(len(avg_ndays)):
        if i == 0:
            avg_ndays[i] == price_arr[0]
        elif i <= n_days:
            avg_ndays[i] = np.average(price_arr[0:i])
        else:
            avg_ndays[i] = np.average(price_arr[i-n_days:i])
    return avg_ndays

# plot closedate and closeprice
def plot_crypto_closes_chart(datetime,price, symbol, plot_avgs = False):
    plt.figure()
    plt.plot(datetime,price)

    if plot_avgs:
        avg50 = get_averages(50,price)
        avg15 = get_averages(15,price)
        avg200 = get_averages(200,price)
        plt.plot(datetime,avg50,label='Avg50')
        plt.plot(datetime,avg15,label='Avg15')
        plt.plot(datetime,avg200,label='Avg200')

    plt.legend([SYMBOL,'Avg50','avg15','avg200'])
    plt.title("Binance Price Chart of "+symbol)
    plt.show()



# ---- START ----

_APIKEY = "wS2uC7zLYtZoB2p7NBznyFek2E3TrsdTNQFtYKj4ESIsNxp2Cm0Pan0tqCzrWw1V"
_SECRET_KEY = "KxqcLBjnqtdxeRxnzFPkIHdV1llULFaHZd0k2Z5k5QgCvGAyawFRnq22VD74Yfms"
_TICKER_API_URL = "https://api.binance.com/api/v1/ticker/price?symbol="
_binsizes = {"1m": 1, "5m": 5, "1h": 60, "1d": 1440}
_batch_size = 750

#test functions
btc_eur_price = get_latest_price("BTCEUR")

#open binance client, some functions
binance_client = Client(api_key=_APIKEY,api_secret=_SECRET_KEY)
depth = binance_client.get_order_book(symbol='BTCEUR')
prices = binance_client.get_all_tickers()
withdraws = binance_client.get_withdraw_history()

# get dataframe of BTCEUR and ETH prices
SYMBOL = "BTCEUR"
SYMBOL2 = 'ETHEUR'
price_table_btc = get_all_binance(SYMBOL,binance_client.KLINE_INTERVAL_1DAY,save=False)
price_table_eth = get_all_binance(SYMBOL2,binance_client.KLINE_INTERVAL_1DAY,save=False)

#extract close date and close price
closes_btc = price_table_btc['close'].to_numpy(dtype=float)
closetime1 = price_table_btc['close_time'].to_numpy()
closedate_btc  = np.array([datetime.utcfromtimestamp(int(ts/1000)) for ts in closetime1])

closes_eth = price_table_eth['close'].to_numpy(dtype=float)
closetime1 = price_table_eth['close_time'].to_numpy()
closedate_eth = np.array([datetime.utcfromtimestamp(int(ts/1000)) for ts in closetime1])



# ------------------ trading szenario--------------------------------------------------
_DAYS_PAST = 100                # start .. days in past
_BUY_FACTOR : float = 0.30               # always buy with 25 % of depot_eur
_SELLFACTOR : float = 0.20              # always sell 25 % of remaining btc

depot_eur : float = 1000                # EUR in depot
depot_btc : float = 0                   # BTC  in depot

days_avg = 7                    # days to build avg of
day_ctr = 0                     # needed in for loop

day_transaction_buy = [0]        # buy dates
price_transaction_buy = [0]      #

day_transaction_sell = [0]       # sell dates
price_transaction_sell = [0]     #

days_buy_sell_pause = 2
day_last_trade = 0
tradeprice_last = 0

depot_value = [0]

prices1 = 30000+2500*np.sin(0.3*np.arange(_DAYS_PAST))

#closes_btc1 = closes_btc[-_DAYS_PAST::]
closes_btc1 = prices1
for btc_price_today in closes_btc1:
    ref_price = np.average(closes_btc[-_DAYS_PAST - days_avg + day_ctr:-_DAYS_PAST + day_ctr])

    if day_last_trade + days_buy_sell_pause < day_ctr:
        if decideBuy(ref_price,btc_price_today):

            buysum = _BUY_FACTOR*depot_eur

            depot_btc += buysum / btc_price_today
            depot_eur -= buysum

            day_transaction_buy.append(day_ctr)
            price_transaction_buy.append(btc_price_today)
            day_last_trade = day_ctr

        elif decideSell(ref_price,btc_price_today):

            sell_btc = _SELLFACTOR*depot_btc
            sellprice = sell_btc*btc_price_today

            depot_eur += sellprice
            depot_btc -= sell_btc

            price_transaction_sell.append(btc_price_today)
            day_transaction_sell.append(day_ctr)
            day_last_trade = day_ctr

    # new reference price is todays price
    day_ctr += 1
    print(btc_price_today,depot_btc,depot_eur)
    depot_value.append(btc_price_today * depot_btc + depot_eur)

print("END")
print(depot_eur)
print(depot_btc)
print(depot_btc * closes_btc1[-1] + depot_eur)

# plot chart
plt.figure()
plt.plot(range(len(closes_btc1)),closes_btc1)
plt.plot(range(len(depot_value)),depot_value)
plt.scatter(np.array(day_transaction_buy),np.array(price_transaction_buy))
plt.scatter(np.array(day_transaction_sell),np.array(price_transaction_sell))
plt.legend(["BTC price","Depot Value","Buy","Sells"])
plt.show()
