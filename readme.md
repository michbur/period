### About
  
Detect periodicity in in quantitative real-time polymerase data.  

**Authors**: [Andrej-Nikolai Spiess](http://www.dr-spiess.de/), [Stefan Roediger](https://www.b-tu.de/fg-multiparameterdiagnostik/projekte/bildbasierte-assays-imagebased-assays-2013-2018/members), [Michal Burdukiewicz](https://github.com/michbur).  

This application is a part of [pcRuniveRsum](https://michbur.github.io/pcRuniveRsum/).

### Reference

[A. Spiess, S. Roediger, M. Burdukiewicz, et al. “System-specific periodicity in quantitative real-time polymerase chain reaction data questions threshold-based quantitation”. Scientific Reports 6 (2016), p. 38951. DOI: 10.1038/srep38951](http://www.nature.com/srep/2016/161213/srep38951/full/srep38951.html).

### Data format

The expected input data is <b>.csv</b> spreadsheet with <b>two</b> columns (see below):

<table border="1" style="width:100%">
  <tr>
    <td><b>Well</b></td>
    <td><b>CQ</b></td> 
  </tr>
  <tr>
    <td>A1</td>
    <td>17.22</td> 
  </tr>
  <tr>
    <td>A2</td>
    <td>17.31</td> 
  </tr>
  <tr>
    <td>...</td>
    <td>...</td> 
  </tr>
</table>

<br>
If the input data is in one row (column-wise), check "Cq values are in one row" below.
 
