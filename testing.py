import numpy as numpy
import plotly.graph_objects as go

d=numpy.zeros(100)
E=numpy.zeros(100)
E1=numpy.zeros(100)
E2=numpy.zeros(100)
for i in range(100):
    d[i]=i/10
    E[i]=(numpy.exp(-numpy.power(i/10,0.5)))
    E1[i]=numpy.exp(-numpy.power(i/10,1))
    E2[i] = 2*numpy.exp(-numpy.power(i/10,0.5))
print(d)
print(E)

Fig = go.Figure(layout=go.Layout(title="Visualisation", uirevision='camera'))
scatt = go.Scatter(x=d,y=E)
scatt1 = go.Scatter(x=d,y=E1)
scatt2 = go.Scatter(x=d,y=E2)

Fig.add_trace(scatt)
Fig.add_trace(scatt1)
Fig.add_trace(scatt2)

Fig.show()
