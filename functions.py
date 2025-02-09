import numpy as np
from lmfit import Parameters, minimize

def read_fp_excelfile(filename, skip_row, cols):
    import pandas as pd
    df=pd.read_excel(filename,
                    header=None,
                    skiprows=skip_row,
                    usecols= cols,
                    names=['Content','FP'])
    #print(df)
    newContent = []
    for item in df['Content']:
        newContent.append(item.split()[1][1:])
    df['Content'] = newContent
    df['Content'] = pd.to_numeric(df['Content'], downcast='integer', errors='coerce')

    data = df.groupby(['Content']).agg(['mean', 'std']).to_numpy()

    return(data)

### Extracting data from specified files
def setup_exps(exps_filename, num_exps):
    import pandas as pd
    fit_param_dict = {}
    if exps_filename != "Default":
        with open(exps_filename) as file:
                for line in file.readlines():
                    if not line.strip().startswith('#'):
                        fit_params = line.strip().split(',')
                        for i in range(len(fit_params)):
                            fit_params[i] = fit_params[i].strip(' ')
                        fit_param_dict.setdefault("Exp",[]).append(fit_params[0])
                        fit_param_dict.setdefault("points",[]).append(int(fit_params[1]))
                        fit_param_dict.setdefault("start",[]).append(int(fit_params[2]))
                        fit_param_dict.setdefault("[titrant]",[]).append(float(fit_params[3]))
                        fit_param_dict.setdefault("dilution",[]).append(float(fit_params[4]))
                        fit_param_dict.setdefault("[fluoro]",[]).append(float(fit_params[5]))
                        fit_param_dict.setdefault("type",[]).append(fit_params[6])

                        if fit_params[6] == 'inhib':
                            fit_param_dict.setdefault("Kd1",[]).append(float(fit_params[7]))
                            fit_param_dict.setdefault("Mt",[]).append(float(fit_params[8]))
                            fit_param_dict.setdefault("units",[]).append(fit_params[9])
                        else:
                            fit_param_dict.setdefault("Kd1",[]).append(0.0)
                            fit_param_dict.setdefault("Mt",[]).append(0.0)
                            fit_param_dict.setdefault("units",[]).append(fit_params[7])

    if exps_filename == "Default":
        for i in range(num_exps):
            fit_param_dict.setdefault("Exp",[]).append("Experiment #" + str(i+1))
            fit_param_dict.setdefault("points",[]).append(int(16))
            fit_param_dict.setdefault("start",[]).append(int(i*16 + 1))
            fit_param_dict.setdefault("[titrant]",[]).append(float(100))
            fit_param_dict.setdefault("dilution",[]).append(float(2.0))
            fit_param_dict.setdefault("[fluoro]",[]).append(float(0.04))
            fit_param_dict.setdefault("type",[]).append("bind")
            fit_param_dict.setdefault("Kd1",[]).append(float(0.0))
            fit_param_dict.setdefault("Mt",[]).append(float(0.0))
            fit_param_dict.setdefault("units",[]).append("μM")

    return pd.DataFrame(fit_param_dict)


def extract_data(data_filename, skip_row, cols):
    if data_filename[-4:] == 'xlsx':
        avg_data = read_fp_excelfile(data_filename, skip_row, cols)
    else:
        exit('Data file should be ".xlsx"')
    return avg_data

def combine_data(avg_data, fit_param_df):
    data_combined = []

    for n in range(len(fit_param_df)):
        name = fit_param_df["Exp"].iloc[n]
        num_points = int(fit_param_df["points"].iloc[n])
        index_start = int(fit_param_df["start"].iloc[n])
        titrant_conc = float(fit_param_df["[titrant]"].iloc[n])
        dil_factor = int(fit_param_df["dilution"].iloc[n])
        fluoro_conc = float(fit_param_df["[fluoro]"].iloc[n])
        type = str(fit_param_df["type"].iloc[n])
        Kd1 = float(fit_param_df["Kd1"].iloc[n])
        Mt = float(fit_param_df["Mt"].iloc[n])
        unit = str(fit_param_df["units"].iloc[n])
    
        y_mean, y_std, x, pt, nKd1, nMt = [], [], [], [], [], []

        conc = titrant_conc
        start = index_start -1
        end = start + num_points
        for i in range(start,end):
            x.extend([conc])
            nKd1.extend([Kd1])
            nMt.extend([Mt])
            y_mean.append(avg_data[i,0])
            y_std.append(avg_data[i,1])
            conc = conc/dil_factor
            pt.append(fluoro_conc)
        
        data_combined.append([name , x , y_mean, y_std, pt, nKd1, nMt, type, unit])

    return data_combined

def setParams(data, concLigand):
# Create fitting parameters for each data set
    fit_params = Parameters()

    fit_params.add('dmax', value=max(data), min=-0.0, max=max(data)*100)
    fit_params.add('dmin', value=min(data), min=-0.0, max=max(data))
    fit_params.add('Kd', value= ((max(concLigand) + min(concLigand))/2), 
                                             min=min(concLigand)/200, max=max(concLigand)*200)

    return(fit_params)

def setParams_inhib(data, concLigand):
# Create fitting parameters for each data set
    fit_params = Parameters()

    fit_params.add('dmax', value=max(data), min=-0.0, max=max(data)*100)
    fit_params.add('dmin', value=min(data), min=-0.0, max=max(data))
    fit_params.add('Kd', value= ((max(concLigand) + min(concLigand))/2), 
                                             min=min(concLigand)/200, max=max(concLigand)*200)

    return(fit_params)

### Functions for fitting a standard binding isotherm
def bindIsotherm(x, dmax, dmin, kd, pt):
    # FP based on Kd for a range of ligand concentrations at a fixed
    # fluorescently labeled protein concentration
    FP = dmin + ( (dmax - dmin) * (((kd + x + pt) - \
                                   (np.sqrt((kd + x + pt)**2 - \
                                            (4 * x * pt))))/(2*pt)))
    return FP

# returns calculated FP at experimental protein and ligand concentrations with fit parameters

def bindIso_dataset(params, concLigand, concProtein):
    """Calculate chemical shift from parameters for data set."""
    dmax = params['dmax']
    dmin = params['dmin']
    Kd = params['Kd']
    return bindIsotherm(concLigand, dmax, dmin, Kd, concProtein)
# returns residuals for optimizing a parameter set for bindIsotherm

def bindIso_objective(params, data, concLigand, concProtein, uncertainty,
                      weight):
    

    resid = (data - bindIso_dataset(params, concLigand, concProtein))

    weighted = np.sqrt(resid ** 2 / uncertainty ** 2)

    if weight:
        return weighted
    else:
        return resid
    
    ### Functions for fitting a competitive inhibitor model
def FPinhib(x, Yml1, Yl1, Kd1, L1t, Kd2, Mt):
    # returns fraction bound for FP inhibition assay
    # Yml1 = dmax, Yl1 = dmin, Kd1 = fixed binding constant, L1t = concProtein (labeled, fixed) 
    # Kd2 = Kd, Mt = protein concentration (fixed, unlabeled)
    L2t = x
    a = (Kd1 + Kd2 + L1t + L2t - Mt)
    b = Kd1*Kd2 + Kd1*(L2t-Mt) + Kd2*(L1t-Mt)
    c = -Mt*Kd1*Kd2
    
    theta = np.arccos( (-2 * a**3 + 9*a*b - 27*c) / (2*np.sqrt((a**2 - 3*b)**3)))
    
    d = ((2 * np.sqrt(a**2 - 3*b) * np.cos(theta/3)) - a)

    FP = Yl1 + ((Yml1-Yl1)*(d/(3*Kd1 + d))) 
    return FP

# returns calculated FP at experimental protein and ligand concentrations with fit parameters
def compInhib_dataset(params, concLigand, concProtein, inhib_Kd1, inhib_Mt):
    """Calculate chemical shift from parameters for data set."""
    dmax = params['dmax']
    dmin = params['dmin']
    Kd = params['Kd']
    return FPinhib(concLigand, dmax, dmin, inhib_Kd1, concProtein, Kd, inhib_Mt)
# returns residuals for optimizing a parameter set for bindIsotherm

def compInhib_objective(params, data, concLigand, concProtein, uncertainty, inhib_Kd1, inhib_Mt,
                      weight):
    

    resid = (data - compInhib_dataset(params, concLigand, concProtein, inhib_Kd1, inhib_Mt))

    weighted = np.sqrt(resid ** 2 / uncertainty ** 2)

    if weight:
        return weighted
    else:
        return resid
    
    ### Fitting data to specified model
def fit_data(data, weight = True):
    from functions import setParams, minimize, bindIso_objective, compInhib_objective
    import numpy as np
    out_params = []

    for i in range(len(data)):
        fit_params = setParams(np.array(data[i][2]), np.array(data[i][1]))

        if data[i][7] == 'bind':
            out = minimize(bindIso_objective, fit_params, method='leastsq', \
                            args=(np.array(data[i][2]), np.array(data[i][1]), \
                                  np.array(data[i][4], dtype='float'), \
                                  np.array(data[i][3], dtype='float'),
                                  weight))
        
        if data[i][7] == 'inhib':
            out = minimize(compInhib_objective, fit_params, method='leastsq', \
                            args=(np.array(data[i][2]), np.array(data[i][1]), \
                                  np.array(data[i][4], dtype='float'), \
                                  np.array(data[i][3], dtype='float'),
                                  np.array(data[i][5], dtype='float'),
                                  np.array(data[i][6], dtype='float'),
                                  weight))

        out_dict = {
                'Kd':out.params['Kd'].value,
                'Kd_err':out.params['Kd'].stderr,
                'dmax':out.params['dmax'].value,
                'dmax_err':out.params['dmax'].stderr,
                'dmin':out.params['dmin'].value,
                'dmin_err':out.params['dmin'].stderr
                }

        out_params.append(out_dict)
    
    return out_params


### Function for plotting data with fit
def pltly(data, out_params,
                num_of_col=2,
                logX = False,
                bgcol = "aliceblue",
                savefile = False,
                savename = 'output',
                savescale = 4,
                writehtml = False,
                htmlfilename = 'output'):
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots

    if len(data) <= num_of_col:
        nrows, ncols = 1, len(data)

        fig = make_subplots(rows=nrows, cols=ncols,
                    shared_xaxes=False,
                    x_title='[Ligand]',
                    y_title='Polarization',
                    )
    else:
        if len(data) % num_of_col == 0:
            nrows, ncols = (len(data)//num_of_col), num_of_col
        else:
            nrows, ncols = ((len(data)//num_of_col)+1), num_of_col
        fig = make_subplots(rows= nrows, cols=ncols,
                    shared_xaxes=False,
                    #vertical_spacing=0.3,
                    x_title='[Ligand]',
                    y_title='Polarization',
                    )

    for i in range(len(data)):
        name, x, y, y_err, pt, nKd1, nMt, type, unit = \
            data[i][0], data[i][1], data[i][2], data[i][3], \
                data[i][4], data[i][5], data[i][6], data[i][7], data[i][8]
        
        Kd, Kd_err, dmax, dmin = \
            out_params[i]['Kd'], out_params[i]['Kd_err'], \
                out_params[i]['dmax'], out_params[i]['dmin']

        row = (i//num_of_col) + 1
        col = (i % num_of_col) + 1 
        
        if logX:
            x_fit = np.logspace(np.log10((min(x))*0.5), np.log10((max(x))*2), 500)
        else:
            x_fit = np.linspace((min(x))*0.9, (max(x))*1.1, 500)

        if type == 'bind':
            pt_fit = np.full(500, pt[0], dtype='float')
            y_fit = bindIsotherm(x_fit, dmax, dmin, Kd, pt_fit)
        elif type == 'inhib':
            Kd1, Mt, L1t = np.full(500, nKd1[0], dtype='float'), np.full(500, nMt[0], dtype='float'), np.full(500, pt[0], dtype='float')
            y_fit = FPinhib(x_fit, dmax, dmin, Kd1, L1t, Kd, Mt)
        
        #fig = go.Figure()

        fig.add_trace(go.Scatter(x = x,y = y, error_y_array = y_err, showlegend = False,
                                     mode = "markers"),
                                     row = row, col = col)
        
        fig.add_trace(go.Scatter(x = x_fit,y = y_fit, showlegend = False,
                                     mode = "lines"),
                                     row = row, col = col)
        fig.add_annotation(text='<b>{0}</b> <br> Kd = {1:.3G} ± {2:.3G} {3}'.format(name,Kd,Kd_err,unit),
                       font=dict(family='Helvetica', size=14),
                       align='center',
                       showarrow=False,
                       xref='x domain',
                       yref='y domain',
                       x=0.5,
                       xanchor='center',
                       y=1.01,
                       yanchor='bottom',
                       row = row, col = col)
        if logX:
            fig.update_xaxes(type='log', row = row, col = col)
    # Update graph style and label axes
    fig.update_xaxes(gridcolor='light gray',gridwidth = 0.2,showgrid = True,title_font_size = 16, mirror = True)
    fig.update_yaxes(gridcolor='light gray',gridwidth = 0.2,showgrid = True,title_font_size = 16, mirror = True)
    fig.update_layout(
                        template = 'simple_white',
                        height = 400 * nrows,
                        width = 500 * ncols,
                        #xaxis = dict(title_text = "[Ligand]"),
                        #yaxis = dict(title_text = "Polarization"),
                        paper_bgcolor=bgcol
                        )
    
    if savefile:  
        import plotly.io as pio
        pio.write_image(fig, savename, scale = savescale)
    
    if writehtml:
        fig.write_html(htmlfilename)
    #fig.show()
    return fig

### Function for plotting data without fit
def pltly_nofit(data,
                num_of_col=2,
                logX = False,
                bgcol = "aliceblue",
                savefile = False,
                savename = 'output',
                savescale = 4,
                writehtml = False,
                htmlfilename = 'output'):
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    
    
    if len(data) <= num_of_col:
        nrows, ncols = 1, len(data)

        fig = make_subplots(rows=nrows, cols=ncols,
                    shared_xaxes=False,
                    x_title='[Ligand]',
                    y_title='Polarization',
                    )
    else:
        if len(data) % num_of_col == 0:
            nrows, ncols = (len(data)//num_of_col), num_of_col
        else:
            nrows, ncols = ((len(data)//num_of_col)+1), num_of_col
        fig = make_subplots(rows= nrows, cols=ncols,
                    shared_xaxes=False,
                    #vertical_spacing=0.3,
                    x_title='[Ligand]',
                    y_title='Polarization',
                    )

    for i in range(len(data)):
        name, x, y, y_err = \
            data[i][0], data[i][1], data[i][2], data[i][3]
        row = (i//num_of_col) + 1
        col = (i % num_of_col) + 1 
        

        fig.add_trace(go.Scatter(x = x,y = y, error_y_array = y_err, showlegend = False,
                                     mode = "markers"),
                                     row = row, col = col)
        
        fig.add_annotation(text='<b>{0}</b>'.format(name),
                       font=dict(family='Helvetica', size=14),
                       align='center',
                       showarrow=False,
                       xref='x domain',
                       yref='y domain',
                       x=0.5,
                       xanchor='center',
                       y=1.01,
                       yanchor='bottom',
                       row = row, col = col)
        if logX:
            fig.update_xaxes(type='log', row = row, col = col)
    # Update graph style and label axes
    fig.update_xaxes(gridcolor='light gray',gridwidth = 0.2,showgrid = True,title_font_size = 16, mirror = True)
    fig.update_yaxes(gridcolor='light gray',gridwidth = 0.2,showgrid = True,title_font_size = 16, mirror = True)
    fig.update_layout(
                        template = 'simple_white',
                        height = 400 * nrows,
                        width = 500 * ncols,
                        #xaxis = dict(title_text = "[Ligand]"),
                        #yaxis = dict(title_text = "Polarization"),
                        paper_bgcolor=bgcol
                        )
    
    if savefile:  
        import plotly.io as pio
        pio.write_image(fig, savename, scale = savescale)
    
    if writehtml:
        fig.write_html(htmlfilename)
    #fig.show()
    return fig

def pltly_table(data,out_params):
    names, Kd, Kd_err, dmax, dmax_err, dmin, dmin_err = [], [], [], [], [], [], []
    for i in range(len(data)):
        names.append(data[i][0])
        Kd.append('{:.2f} ± {:.2f}'.format(out_params[i]['Kd'], out_params[i]['Kd_err']))
        dmax.append('{:.2f} ± {:.2f}'.format(out_params[i]['dmax'], out_params[i]['dmax_err']))
        dmin.append('{:.2f} ± {:.2f}'.format(out_params[i]['dmin'], out_params[i]['dmin_err']))

    import plotly.graph_objects as go
    headerColor = 'grey'
    rowEvenColor = 'lightgrey'
    rowOddColor = 'white'
    
    fig = go.Figure(data=[go.Table(
        columnwidth = [200,100],
        header=dict(values=['Exp', 'Kd', 'dmax', 'dmin'],
                    line_color='darkslategray',
                    fill_color=headerColor,
                    align=['left','center'], # 'center',
                    font=dict(color='white', size=12)
                    ),
        cells=dict(values=[ names, Kd, dmax, dmin],
                    # line_color='darkslategray',
                    line_color= [[rowOddColor,rowEvenColor,rowOddColor, rowEvenColor,rowOddColor]*5],
                    #fill_color='white',
                    fill_color = [[rowOddColor,rowEvenColor,rowOddColor, rowEvenColor,rowOddColor]*5],
                    align=['left','center'], # 'center',
                    font = dict(color = 'darkslategray', size = 11),
                    )
                    )
    ])

    fig.update_layout(width=800, height=400)
    #fig.show()
    return fig

def create_table(data,out_params):
    import pandas as pd
    table_dict = {}
    ind = []
    for i in range(len(data)):
        table_dict.setdefault("Exp",[]).append(data[i][0])
        table_dict.setdefault("Kd",[]).append('{:.2f} ± {:.2f} {}'.format(out_params[i]['Kd'], out_params[i]['Kd_err'], data[i][8]))
        table_dict.setdefault("FPmax",[]).append('{:.2f} ± {:.2f}'.format(out_params[i]['dmax'], out_params[i]['dmax_err']))
        table_dict.setdefault("FPmin",[]).append('{:.2f} ± {:.2f}'.format(out_params[i]['dmin'], out_params[i]['dmin_err']))
    return pd.DataFrame(table_dict)

def make_tables(data):
    import pandas as pd
    import panel as pn
    from bokeh.models.widgets.tables import NumberFormatter
    pn.extension('tabulator')

    num_formatter = {
        '[ligand]' : NumberFormatter(format='0.00'),
        'FP' : NumberFormatter(format='0.00'),
        'FP error' : NumberFormatter(format='0.00')
    }
    tables = pn.Row()
    tables_df =[]
    for n in range(len(data)):
        table_df = pd.DataFrame({
            '[ligand]': data[n][1],
            'FP': data[n][2],
            'FP error': data[n][3]})
        table = pn.widgets.Tabulator(table_df, 
                                    sortable=False,
                                    #show_index = False,
                                    theme = 'simple',
                                    layout = 'fit_data',
                                    header_align = 'center',
                                    disabled = True,
                                    show_index = False,
                                    formatters = num_formatter,
                                    text_align = 'right'
                                    )
        title_text = pn.widgets.StaticText(name='Table', value = data[n][0])
        col = pn.Column(title_text,table)
        tables.append(col)
        tables_df.append(table_df)
    return tables, tables_df