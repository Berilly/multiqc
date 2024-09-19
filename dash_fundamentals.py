# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

from dash import Dash, html, dcc
import plotly.express as px
import pandas as pd

multiqc_df = pd.read_csv("QC-2024-09-19-NextSeq_multiQC.csv", encoding="latin1")
multiqc_df_failed = (multiqc_df[multiqc_df["Yield"] == "FAILED"])
statistics_df = pd.read_csv("QC-2024-09-19-NextSeq_stat.csv", encoding="latin1")
print(statistics_df)

def generate_table(dataframe, max_rows=10):
    return html.Table([
        html.Thead(
            html.Tr([html.Th(col) for col in dataframe.columns])
        ),
        html.Tbody([
            html.Tr([
                html.Td(str(dataframe.iloc[i][col])) for col in dataframe.columns
            ]) for i in range(min(len(dataframe), max_rows))
        ])
    ])

app = Dash(__name__)

colors = {
    'background': '#FFA500',
    'text': '#7FDBFF'
}

fig1 = px.sunburst(multiqc_df,
                  path=['Yield', 'ProjectName'], 
                  values='X..Reads',
                  title="Sunburst Plot for All Data",
                  color='ProjectName', 
                  color_continuous_scale='rainbow',
                  )  
fig1 = fig1.update_traces(maxdepth=3)

fig2 = px.sunburst(multiqc_df_failed,
                  path=['Yield', 'ProjectName', 'Isolation Kit', 'SampleID'], 
                  values='X..Reads',
                  title="Sunburst Plot for Failed Data",
                  color='ProjectName', 
                  color_continuous_scale='rainbow',
                  )  
fig2 = fig2.update_traces(maxdepth=3)

app.layout = html.Div(style={'backgroundColor': colors['background']}, children=[
    html.H1(
        children='Sequrencing Data Analysis',
        style={
            'textAlign': 'center',
            'color': colors['text'],
        }
    ),
    html.Div(
        children='Results of the last run', 
        style={
            'textAlign': 'center',
            'color': colors['text']
        }
    ),
    html.Div(
    html.H4(children = generate_table(statistics_df)
            )
    ),
    dcc.Graph(
        id='example-graph',
        figure=fig1
    ),
    dcc.Graph(
        id='example-graph2',
        figure=fig2
    )
])

if __name__ == '__main__':
    app.run(debug=True)
