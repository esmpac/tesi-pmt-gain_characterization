import ROOT                 
from IPython.display import Image, display

# Impedisce a ROOT di aprire finestre GUI (utile in Colab)
ROOT.gROOT.SetBatch(True)

from pathlib import Path
import numpy as np
import pandas as pd
import time


VOLTAGE_MIN = 900
VOLTAGE_MAX = 1500


CALIBRAZIONE_ADC=0.0132
CARICA_ELETTRONE=1.6*10**-7 #in pC


ROOT.gStyle.SetOptFit(1111)
CHANNELS = range(7)

pedestal= [882.91, 886.25, 887.94, 882.36, 880.55, 887.07, 883.53]

def ped(x, par):                                                
    gauss0 = par[0]*np.exp(-0.5*((x[0])/par[1])**2)
    return float(gauss0)




    
def double_gaus(x, par):
    #return float(par[0]*np.exp(-0.5*((x[0]-par[1])/par[2])**2))
    #exp = par[0]*np.exp(-1*x[0]/par[1])
    gauss0 = par[0]*np.exp(-0.5*((x[0])/par[1])**2)
    gauss1 = par[2]*np.exp(-0.5*((x[0]-par[3])/par[4])**2)
    #gauss2 = par[5]*np.exp(-0.5*((x[0]-2*par[3])/par[4])**2)
    #return float(gauss0+gauss1+gauss2)
    return float(gauss0+gauss1)


def PEDESTAL():               

    info_pedestal = {}
    hist_list = [] 

    pedestal_data = Path("/content/drive/MyDrive/tesi_swgo_2025/tesi-pmt-gain_characterization/data/batch_3/pedestal_characterization/2025_04_09/acq_1")   
    print("Directory pedestal:", pedestal_data)                     

    # Prepariamo il fit gaussiano
    fit_ped = ROOT.TF1("f", "gaus", 875, 895, 3)                      
    fit_ped.SetParNames("A", "mu", "sigma")

    # Lista per accumulare tutti i CSV
    df_list = []
    for data_file in pedestal_data.rglob("*.csv"):
        print("Caricamento file:", data_file)
        df_list.append(pd.read_csv(data_file))

    # Verifica che ci siano file
    if not df_list:
        raise FileNotFoundError(f"Nessun file CSV trovato in {pedestal_data}")

    # Concatenazione di tutti i CSV in un unico DataFrame
    df = pd.concat(df_list, ignore_index=True)

    # Canvas ROOT per il plot
    canvas = ROOT.TCanvas("a", "Multi-Channel Histogram Pedestal", 1800, 1000)          
    canvas.Divide(3, 3)                                                        

    for channel in CHANNELS:                                                        

        hist = ROOT.TH1D(f"hist_ch{channel}_pedestal", f"Channel {channel}", 30, 870, 900)    
        hist_list.append(hist)     

        energy_values = df[df["Channel"] == channel]["Energy"].values

        for ener in energy_values:                                                   
            hist.Fill(ener)                                                         
  
        print("Pedestal Channel:", channel)

        canvas.cd(channel+1)
        hist.SetLineColor(ROOT.kBlue)
        hist.SetFillColorAlpha(ROOT.kBlue, 0.3)

        entries = hist.GetEntries()
        mean    = hist.GetMean()
        std_dev = hist.GetStdDev()
        k = 4
    
        fit_ped.SetParameters(int(entries/100), mean, std_dev)
        fit_ped.SetRange(mean - k*std_dev, mean + k*std_dev)

        hist.Fit(fit_ped, "RL")
        n_ped = fit_ped.GetParameter(0)
        mu_ped = fit_ped.GetParameter(1)
        sigma_ped = fit_ped.GetParameter(2)

        print(f"A: {n_ped}, mu: {mu_ped}, sigma: {sigma_ped}")
    
        hist.Draw() 
        info_pedestal[channel] = mu_ped
        
            # Mostra grafico in Colab
        canvas_channel = ROOT.TCanvas(f"c_ped_{channel}", f"Channel {channel}", 800, 600)
        hist.Draw()
        fit_ped.Draw("same")
        filename = f"pedestal_ch{channel}.png"
        canvas_channel.SaveAs(filename)
        display(Image(filename))

    canvas.Modified() 
    canvas.Update()

    print(f"I valori medi dei piedistalli: {info_pedestal}")
    print("Press any key to stop the program...")
    input()

    return info_pedestal                           

'''
def PEDESTAL():               


    info_pedestal = {}
    hist_list = [] 

    pedestal_data = Path("/content/drive/MyDrive/tesi_swgo_2025/tesi-pmt-gain_characterization/data/batch_3/pedestal_characterization/2025_04_09/acq_1")   
    print(pedestal_data)                     

    fit_ped = ROOT.TF1("f", "gaus", 875, 895, 3)                      
    fit_ped.SetParNames("A", "mu", "sigma")


    for data_file in pedestal_data.rglob("*.csv"): 
        print(data_file)               
        df = pd.read_csv(data_file)                                
            

    canvas = ROOT.TCanvas("a", "Multi-Channel Histogram Pedestal", 1800, 1000)          
    canvas.Divide(3, 3)                                                        


    for channel in CHANNELS:                                                        

        hist = ROOT.TH1D(f"hist_ch{channel}_pedestal", f"Channel {channel}", 30, 870, 900)    

        hist_list.append(hist)     

        energy_values = df[df["Channel"] == channel]["Energy"].values


        for ener in energy_values:                                                   
            hist.Fill(ener)                                                         
  
        print("Pedestal Channel:", channel)

    
        canvas.cd(channel+1)
        hist.SetLineColor(ROOT.kBlue)
        hist.SetFillColorAlpha(ROOT.kBlue, 0.3)


        entries = hist.GetEntries()
        mean    = hist.GetMean()
        std_dev = hist.GetStdDev()
        k = 4
    
    
        fit_ped.SetParameters(int(entries/100), mean, std_dev)
        fit_ped.SetRange(mean - k*std_dev, mean + k*std_dev)

        hist.Fit(fit_ped, "RL")
        n_ped = fit_ped.GetParameter(0)
        mu_ped = fit_ped.GetParameter(1)
        sigma_ped = fit_ped.GetParameter(2)

        print("A:", n_ped, "\n mu:", mu_ped, "\n sigma:",sigma_ped)
    
                                    
        hist.Draw() 

        info_pedestal[channel] = mu_ped                                                                  


    canvas.Modified() 
    canvas.Update()

    print(f"I valori medi dei piedistalli: {info_pedestal}")

    print("Press any key to stop the programm...")
    input()

    return info_pedestal
'''




def SPE(pedestal_values):

    
    hist_list = []
    guadagno_flt = []
    guadagno = []


    fit_ped = ROOT.TF1("f", ped, 25, 80, 2)
    fit_ped.SetParNames("A", "sigma0")

    fit_prima = ROOT.TF1("s", "gaus", 100, 500, 3)
    fit_prima.SetParNames("B, mu, sigma")



    spe_data_long = Path("/content/drive/MyDrive/tesi_swgo_2025/tesi-pmt-gain_characterization/data/batch_3/single_photoelectron/2025_04_10")


    for data_file in spe_data_long.rglob("*.csv"):
        df = pd.read_csv(data_file)

    
    canvas = ROOT.TCanvas("a", "Multi-Channel Histogram SPE", 1800, 1000)
    canvas.Divide(3, 3)

    
    for channel in CHANNELS:

        hist = ROOT.TH1D(f"hist_ch{channel}_spe", f"Channel {channel}", 500, 0, 1000)
        hist_list.append(hist)


        energy_values = df[df["Channel"] == channel]["Energy"].values


        for ener in energy_values:
            hist.Fill(ener-pedestal_values[channel])      


        
        canvas.cd(channel+1)

        hist.SetLineColor(ROOT.kBlue)
        hist.SetFillColorAlpha(ROOT.kBlue, 0.3)

        max_scale = int(hist.GetEntries()/10)

        ##################################
        fit_ped.SetParameters(max_scale, 40)
        fit_ped.SetParLimits(0, 0, max_scale) #A ped
        fit_ped.SetParLimits(1, 0, 60) #sigma ped
        hist.Fit(fit_ped, "RL")

        #################################

        ##############################
        n_ped = fit_ped.GetParameter(0)
        sigma_ped = fit_ped.GetParameter(1)
        ####################################





        fit_prima.SetParameters(max_scale, 300, 40)


        hist.Fit(fit_prima, "RL")

        n_prima = fit_prima.GetParameter(0)
        mu_prima = fit_prima.GetParameter(1)
        sigma_prima = fit_prima.GetParameter(2)

        if mu_prima <= 200:
            fit_function = ROOT.TF1("f", double_gaus, 35, 300, 5)
            fit_function.SetParNames("A", "sigma0", "B", "mu", "sigma")

        
            opt = ROOT.Fit.DataOptions()
            rangeB = ROOT.Fit.DataRange()
            rangeB.SetRange(35, 300)
        else:
            fit_function = ROOT.TF1("f", double_gaus, 35, 500, 5)
            fit_function.SetParNames("A", "sigma0", "B", "mu", "sigma")

        

        fit_function.SetParameters(n_ped, sigma_ped, n_prima, mu_prima, sigma_prima)

        
        fit_function.SetParLimits(2, 0.5*n_prima, 1.5*n_prima) #B prima gauss
        fit_function.SetParLimits(3, 0.5*mu_prima, 1.5*mu_prima) #mu 
        fit_function.SetParLimits(4, 0.5*sigma_prima, 1.5*sigma_prima) #sigma


        
        hist.Fit(fit_function, "R")

        guadagno_flt.append((fit_function.GetParameter(3) * CALIBRAZIONE_ADC))

        valore = format((fit_function.GetParameter(3) * CALIBRAZIONE_ADC)/ CARICA_ELETTRONE, ".2e")
        guadagno.append(valore)
        
        # ---------------------- NUOVO ----------------------
        canvas_channel = ROOT.TCanvas(f"c_spe_{channel}", f"SPE Channel {channel}", 800, 600)
        hist.Draw()
        fit_function.Draw("same")
        filename = f"spe_ch{channel}.png"
        canvas_channel.SaveAs(filename)
        display(Image(filename))

    
    print(guadagno)  
    canvas.Modified() 
    canvas.Update()
    print("Press any key to stop the programm...")
    input()

    return guadagno_flt






'''

def GUADAGNO(pedestal_values):
    info_gaudagno = {}


    guadagno_data = Path("/content/drive/MyDrive/tesi_swgo_2025/tesi-pmt-gain_characterization/data/batch_3/gain_curve/2025_04_10")

    for data_file in guadagno_data.rglob("*.csv"):

        try:
            voltage = int(data_file.name.split("_")[-1].split(".")[0])
        except:
            print("Impossibile recuperare il valore della tensione")
            return

        canvas = ROOT.TCanvas(f"{voltage}", f"Multi-Channel Histogram- V {voltage}", 1800, 1000)
        canvas.Divide(3, 3)

        hist_list = []
        print(voltage)
        
        df = pd.read_csv(data_file)

        for channel in CHANNELS:

            hist = ROOT.TH1D(f"hist_ch{channel}_v_{voltage}", f"Channel {channel} ", 1000, 0, 16000)
            hist_list.append(hist)

            channel_data = df[df["Channel"] == channel]["Energy"]

            if not channel_data.empty:
                print(channel_data)

                for ener in channel_data.values:
                    hist.Fill(ener-pedestal_values[channel])  

                
                mean = hist.GetMean()
                sigma = hist.GetRMS()
                amplitude = int(hist.GetEntries()/10)
                min_value = mean - 5*sigma
                max_value = mean + 5*sigma
                hist.SetAxisRange(min_value, max_value, "X")

                fit_gauss = ROOT.TF1("f", "gaus", min_value, max_value, 3)
                fit_gauss.SetParameters(amplitude, mean, sigma)

                if voltage not in info_gaudagno:
                    info_gaudagno[voltage] = []  

                info_gaudagno[voltage].append((channel, fit_gauss.GetParameter(1)* CALIBRAZIONE_ADC))
                
                

            canvas.cd(channel+1)
            hist.SetLineColor(ROOT.kBlue)
            hist.SetFillColorAlpha(ROOT.kBlue, 0.3)
            hist.Fit(fit_gauss, "R")
            # ---------------------- NUOVO ----------------------
            canvas_channel = ROOT.TCanvas(f"c_gain_{voltage}_ch{channel}", f"Gain Ch {channel} V {voltage}", 800, 600)
            hist.Draw()
            fit_gauss.Draw("same")
            filename = f"gain_{voltage}_ch{channel}.png"
            canvas_channel.SaveAs(filename)
            display(Image(filename))
            #hist.Draw("")

        canvas.Modified() 
        canvas.Update()
        print("Press any key to stop the programm...")
        input()
    
    return info_gaudagno
'''

def GUADAGNO(pedestal_values):
    info_gaudagno = {}

    guadagno_data = Path("/content/drive/MyDrive/tesi_swgo_2025/tesi-pmt-gain_characterization/data/batch_3/gain_curve/2025_04_10")

    for data_file in guadagno_data.rglob("*.csv"):

        try:
            voltage = int(data_file.name.split("_")[-1].split(".")[0])
        except:
            print("Impossibile recuperare il valore della tensione")
            return

        canvas = ROOT.TCanvas(f"{voltage}", f"Multi-Channel Histogram- V {voltage}", 1800, 1000)
        canvas.Divide(3, 3)

        hist_list = []
        print(voltage)
        
        df = pd.read_csv(data_file)

        for channel in CHANNELS:

            hist = ROOT.TH1D(f"hist_ch{channel}_v_{voltage}", f"Channel {channel} ", 1000, 0, 16000)
            hist_list.append(hist)

            channel_data = df[df["Channel"] == channel]["Energy"]

            if not channel_data.empty:
                print(channel_data)

                for ener in channel_data.values:
                    hist.Fill(ener-pedestal_values[channel])  

                mean = hist.GetMean()
                sigma = hist.GetRMS()
                amplitude = int(hist.GetEntries()/10)
                min_value = mean - 5*sigma
                max_value = mean + 5*sigma
                hist.SetAxisRange(min_value, max_value, "X")

                fit_gauss = ROOT.TF1("f", "gaus", min_value, max_value, 3)
                fit_gauss.SetParameters(amplitude, mean, sigma)

                if voltage not in info_gaudagno:
                    info_gaudagno[voltage] = []  

                hist.Fit(fit_gauss, "R")
                info_gaudagno[voltage].append((channel, fit_gauss.GetParameter(1)* CALIBRAZIONE_ADC))
                
                # ---------------------- NUOVO ----------------------
                canvas_channel = ROOT.TCanvas(f"c_gain_{voltage}_ch{channel}", f"Gain Ch {channel} V {voltage}", 800, 600)
                hist.SetLineColor(ROOT.kBlue)
                hist.SetFillColorAlpha(ROOT.kBlue, 0.3)
                hist.Draw()
                fit_gauss.Draw("same")
                filename = f"gain_{voltage}_ch{channel}.png"
                canvas_channel.SaveAs(filename)
                display(Image(filename))

        canvas.Modified() 
        canvas.Update()
        print("Press any key to stop the programm...")
        input()
    
    return info_gaudagno





def GAIN_CURVE(pedestal_values):
    #ROOT.gROOT.SetBatch(True)
    carica_spe = SPE(pedestal_values)
    carica_raw = GUADAGNO(pedestal_values)
    #ROOT.gROOT.SetBatch(False)


    npe = []
    tensioni_gaudagno = []

    gain_curve = {}

    carica_1450 = carica_raw[1450]

    for data in carica_1450:
        channel = data[0]
        charge = data[1]
        npe_1450 = charge / (carica_spe[channel])
        npe.append(npe_1450)


    for voltage in range(900, 1500, 50):

        carica = carica_raw[voltage]

        for data in carica:

            channel = data[0]
            charge = data[1]

            gaudagno = charge / (npe[channel] * CARICA_ELETTRONE)

            if channel not in gain_curve:
                gain_curve[channel] = []  

            gain_curve[channel].append((voltage, gaudagno))


    print(npe)
    canvas = ROOT.TCanvas(f"Gain Curves", f"Gain Curves", 1800, 1000)
    canvas.Divide(3, 3)
    

    graph_list = []
    line_list = []
    vertical_lines = []

    for channel in CHANNELS:

        num_values = len(gain_curve[channel])
        graph = ROOT.TGraph(num_values)

        linea = ROOT.TLine(VOLTAGE_MIN, 3.5e6, VOLTAGE_MAX, 3.5e6)
        line_list.append(linea)
        

        for data in gain_curve[channel]:
            index = gain_curve[channel].index(data)
            voltage = data[0]
            gain = data[1]
            graph.SetPoint(index, voltage, gain)
            

            graph.SetTitle(f"Channel: {channel}")
            graph_list.append(graph)

        power = ROOT.TF1(f"power_ch{channel}", "[0]*pow(x, [1])", VOLTAGE_MIN, VOLTAGE_MAX, 2)
        power.SetParameters(1e-14, 4)
        power.SetLineColor(ROOT.kBlack)

        
        
        canvas.cd(channel+1)
        ROOT.gPad.SetLogy(1)
        graph.SetMarkerStyle(20)  
        graph.SetMarkerColor(ROOT.kRed)
        graph.Draw("AP")

        graph.Fit(power, "R")

        intersezione = power.GetX(3.5e6, VOLTAGE_MIN, VOLTAGE_MAX)
        tensioni_gaudagno.append(intersezione)
        print(f"L'intersezione per il canale: {channel} è a x =", intersezione)

        linea.SetLineColor(ROOT.kBlue)  
        linea.SetLineWidth(2)
        linea.Draw("same")
        
        
        
        
        canvas_channel = ROOT.TCanvas(f"c_gain_curve_ch{channel}", f"Gain Curve Channel {channel}", 800, 600)
        ROOT.gPad.SetLogy(1)
        graph.Draw("AP")
        power.Draw("same")
        linea.Draw("same")
        
        # Creo e disegno vertical_line SOLO DOPO averla definita
        vertical_line = ROOT.TLine(intersezione, 1e6, intersezione, 1e7)
        vertical_line.SetLineColor(ROOT.kGreen)
        vertical_line.SetLineWidth(2)
        vertical_line.Draw("same")
        vertical_lines.append(vertical_line)
        
        filename = f"gain_curve_ch{channel}.png"
        canvas_channel.SaveAs(filename)
        display(Image(filename))
        
        
        
        
        
        
        
        
        
        
        
    '''
        canvas_channel = ROOT.TCanvas(f"c_gain_curve_ch{channel}", f"Gain Curve Channel {channel}", 800, 600)
        ROOT.gPad.SetLogy(1)
        graph.Draw("AP")
        power.Draw("same")
        linea.Draw("same")
        vertical_line.Draw("same")
        filename = f"gain_curve_ch{channel}.png"
        canvas_channel.SaveAs(filename)
        display(Image(filename))

        vertical_line = ROOT.TLine(intersezione, 1e6, intersezione, 1e7)
        vertical_lines.append(vertical_line)
        vertical_line.SetLineColor(ROOT.kGreen)
        vertical_line.SetLineWidth(2)
        vertical_line.Draw("same")
   '''
   
    canvas.Modified() 
    canvas.Update()
    print("Press any key to stop the programm...")
    input()





            


if __name__ == '__main__':

    #a = SPE()
    # print(a)

    piedistallo = PEDESTAL()
    GAIN_CURVE(piedistallo)
    #SPE()       
    #GUADAGNO()
    #PEDESTAL()
    #print (pedestal_values_after_fit) #change
    #GAIN_CURVE()
