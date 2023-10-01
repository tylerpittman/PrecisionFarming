# python3 -tt /Users/tylerpittman/Farm/agYieldProject/weatherDownloadScript.py
#https://api.weather.com/v2/pws/history/all?stationId=ILACADEN2&format=json&units=m&date=20180701&apiKey=

import requests
import pandas as pd

station_ID = 'ILACADEN2'
station_key = ''
start_date = 20180701 
end_date = 20190630

for x in range(start_date, end_date + 1):
    WU_date = str(x)
    
    url = 'https://api.weather.com/v2/pws/history/all?stationId=' + station_ID + '&format=json&units=m&date=' + WU_date + '&apiKey=' + station_key
    df = pd.read_json(url)
    tf = pd.read_json( (df['observations']).to_json(), orient='index')
    tf2 = pd.read_json( (tf['metric']).to_json(),  orient='index')
    tf3 = tf.join(tf2)
    
    tf_out = tf3[['obsTimeLocal', 'tempHigh', 'dewptHigh', 'winddirAvg', 'windspeedAvg','windgustHigh','pressureMax',
                  'humidityAvg', 'precipRate', 'precipTotal', 'uvHigh', 'solarRadiationHigh']]
    
    tf_sum = {'Record' :['obsTimeLocal', 'tempHigh', 'tempLow', 'dewptHigh', 'dewptLow', 'winddirAvg', 'windspeedAvg','windgustHigh',
                         'pressureMax', 'pressureMin', 'humidityAvg','precipRate', 'precipTotal', 'uvHigh', 'solarRadiationHigh'],
            'Observation' :[tf3['obsTimeLocal'].max(), tf3['tempHigh'].max(), tf3['tempLow'].min(), tf3['dewptHigh'].max(), 
                            tf3['dewptLow'].min(), tf3['winddirAvg'].mean(), tf3['windspeedAvg'].mean(), tf3['windgustHigh'].max(), 
                            tf3['pressureMax'].max(), tf3['pressureMin'].min(), tf3['humidityAvg'].mean(), tf3['precipRate'].max(), 
                            tf3.at[275, 'precipTotal'], tf3['uvHigh'].mean(), tf3['solarRadiationHigh'].mean()]
    }
    
    tf_sum_base = pd.DataFrame(tf_sum, columns=['Record', 'Observation'])
    tf_sum_out = tf_sum_base.transpose()
    
    print(WU_date)
    
    out_path = "/Users/tylerpittman/Farm/agYieldProject/data" + WU_date + ".xlsx"
    with pd.ExcelWriter(out_path) as writer:
        tf_out.to_excel(writer , sheet_name= WU_date)
        tf_sum_out.to_excel(writer , sheet_name= 'Summary')
