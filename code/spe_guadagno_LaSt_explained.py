import ROOT                                                 # importo il modulo pyROOT che:

                                                            # ponendosi come "ponte" fra il framework ROOT già istallato e scritto in C++ permette di usare tutte 
                                                            # le sue funzionalità dentro Python. per il calcolo scientifico e per l'analisi dati
                                                            
                                                            # sa TRADURRE le chiamate Python in chiamate alle classi C++ di ROOT (chat gpt fondo chat https://chatgpt.com/c/68e005ad-7308-8328-935d-00e989df0562)

                                                            # carica le librerie del framework ROOT solo quando servono ad esempio per 
                                                            # la definizione di funzioni di fit (TF1), la gestione della grafica tramite TCanvas, la creazione di istogrammi (TH1D), ecc....

                                                           


from pathlib import Path                                    # importo dal modulo pathlib, un modulo della standard library di Python per la gestione
                                                            # dei percorsi dei file e delle cartelle in modo più efficiente rispetto alle stringhe (obsolete), la classe Path, che 
                                                            # offre numerosi metodi per la manipolazione dei percorsi in modo semplice e leggibile
                                                            # In questo modo potremo creare oggetti della classe Path e operare su di essi tramite dei metodi della classe Path
                                                           
                                                            # Quello che succede è che si porrà una variabile uguale al percorso del file facendolo quindi divenire un oggetto della classe Path 
                                                            # su cui potremo operare con i metodi della classe Path


import numpy as np                                          # numpy è una libreria per il calcolo numerico in Python,
                                                            # particolarmente efficiente per la gestione di array e operazioni matematiche avanzate.

import pandas as pd                                         # pandas è una libreria per la manipolazione e l'analisi di dati strutturati in Python,
                                                            # utile per leggere e organizzare i dati da file .csv.

import time                                                 # time è una libreria della standard library di Python per gestire operazioni legate al tempo.








VOLTAGE_MIN = 900                                                               # è il voltaggio minimo a cui poniamo i diversi PMT
VOLTAGE_MAX = 1500                                                              # è il voltaggio massimo a cui poniamo i diversi PMT


CALIBRAZIONE_ADC=0.0132                                                         # è la costante di conversione per passare dal canale ADC a cui visulizzo il segnale di interesse alla carica di quel segnale 
                                                                                # ( CONVERTO CANALE ADC ----> CARICA IN pC)
CARICA_ELETTRONE=1.6*10**-7                                                     # carica elettrone

colors = ["green", "magenta", "red", "orange", "purple", "blue", "brown"]       # definisce una lista di colori, rappresentati come stringhe

ROOT.gStyle.SetOptFit(1111)                                                     # serve a visualizzare i risultati del fit direttamente sul grafico in ROOT
                                                                                # SetOptFit accetta un intero il cui valore (bitwise) decide quali informazioni del fit vengono stampate sul canvas (di solito in un box in alto a destra).
                                                                                # Il numero 1111 è una combinazione di bit, ciascuno dei quali attiva una specifica informazione
                                                                                # Quindi SetOptFit(1111) oppure SetOptFit(15) fa apparire:

                                                                                # Nome della funzione
                                                                                # Parametri del fit (con errori)
                                                                                # Chi² / NDF
                                                                                # P-value (Probabilità del fit)


BATCH = 1                                                                       # la batch è posta essere uguale ad 1 ( costante ) anche se vale che la batch è 1 sola

CHANNELS = range(7)                                                             # crea un 'oggetto range' che rappresenta una sequenza di numeri interi da 0 a 6 (cioè 7 numeri in totale, escluso il 7)
                                                                                # tecnicamente range(7) non è una lista, è un oggetto più leggero in termini di memoria, perché non crea subito tutti gli elementi, 
                                                                                # ma li genera al volo quando servono (è un oggetto iterabile). E' come se fosse un qualcosa che è pronto ad essere iterato con un 
                                                                                # ciclo for ma non è un lista vera e propria a meno che non cicli su di essa con un cico for e stampi a video l'output 
                                                                                # r = range(5)         # Oggetto range                  
                                                                                # lista = list(r)      # Ora lista è una vera lista: [0, 1, 2, 3, 4]

                                                                                # print(lista)         # Output: [0, 1, 2, 3, 4]
                                                                                # print(type(lista))   # <class 'list'>

pedestal= [882.91, 886.25, 887.94, 882.36, 880.55, 887.07, 883.53]              # piedistallli che vengono sottratti come off-set al segnale rivelato











def ped(x, par):                                                                   # gaussiana ( in verità è un esponenziale ) di fit / funzione di fit gaussiana per il piedistallo --> rappresenta la forma matematica che vogliamo usare per il fit sul piedistallo
    gauss0 = par[0]*np.exp(-0.5*((x[0])/par[1])**2)                                # ha due parametri ed è CENTRATA IN ZERO ( li inizializzo prima del fit globale sugli istogrammi su SPE )
                                                                                   # è centrata in 0 perchè ai valori in energia per i vari PMT come si vede sotto sottraiamo i valori del pedestal e quindi 
                                                                                   # l'offset in energia che è dovuto al rumore/ noise dell'elettronica, a PMT non alimentato.
                                                                                   # Sottratto l'offset, e quindi il valor medio della distribuzione gaussiana del rumore diversa per ciascun PMT non alimentato, 
                                                                                   # varrà che l'espoenziale deve partire da 0 ( NELLA IPOTESI CHE L'OFFSET IN ENERGIA A PMT NON ALIMENTATO E A PMT ALIMENTATO SIANO COMPATIBILI) e quindi conterrà solo due parametri e non 3.
                                                                                   # In questo senso, per automatizzare il metodo, questo modo di procedere è eccellente, perchè avrò sempre una discesa espenenziale che parte da 0 
                                                                                   # identificata da 2 parametri e poi la curva di SPE.
                                                                                   # Se trovo che   il segnale di 2PE - 2 volte il segnale di SPE è compatibile con lo 0   allora avendo operato la traslazione, avendo sottratto per i valori del pedestal, varrà che 
                                                                                   # automaticamente i valori del pedestal a PMT non alimentato e PMT alimentato sono compatibili ( in questo senso ho solo traslato tutto lungo l'asse x e poi dal fit sugli istogrammi visto che   il segnale di 2PE - 2 volte il segnale di SPE è compatibile con lo 0 
                                                                                   # ho mostrato che l'ipotesi che gli offset siano compatibili è corretta)

                                                                                   # E' importante ribadire che la THRESHOLD IN ENERGIA è tale che vedo solo la parte finale della coda destra del picco gaussiano 
                                                                                   # di rumore che tratto come un esponenziale decrescente che può essere fittato a partire dai PRIMI VALORI in cui " ho effettivamente
                                                                                   # i primi punti da fittare " nel caso in cui sottragga il valore del pedestal.
                                                                                   # E quindi quello che accade è che per ogni PMT sto solo sottraendo l'offeset relativo (pedestal) dai dati ricavati per poi fare un 
                                                                                   # fit sopra essi --> è una traslazione lungo l'asse x
                                                                                   # il MODELLO che fitto dopo aver sottratto l'offset è quello di un esponenziale decrescente che parte da 0 ma il range in cui andrò a
                                                                                   # fittare il modello inizia da dove ho i punti per poter fare il fit, nonostante il modello parta da 0.
                                                                                   # (lavoro nello " spazio dei dati a cui ho sottratto il pedestal ")
    #gauss0 = par[0]*np.exp(-1*x[0]/par[1])
    return float(gauss0)                                                           # il return float(gauss0) è necessario perchè la funzione 'ped' ritorni un float e quindi serve a garantire che il 
                                                                                   # numero restituito dalla funzione ped sia un numero di tipo 'float'. In Python, quando calcoli un'espressione come 
                                                                                   # np.exp(-0.5*((x[0])/par[1])**2), il risultato potrebbe essere di tipo numpy.float64 o simile, che è un tipo 
                                                                                   # numerico ma non esattamente un float standard di Python (che è di tipo float in Python, ovvero float64 in numpy).
                                                                                   # Quindi con quella riga forza la conversione esplicita del risultato in un numero reale Python

                                                                                   # La conversione esplicita a float (con la funzione float()) serve a evitare problemi di compatibilità o di tipo in 
                                                                                   # altre parti del codice che potrebbero aspettarsi un float nativo di Python e non un tipo di dato numpy.




def primm(x, par):                                                                 # gaussiana di fit / funzione di fit gaussiana per il single-photo-electron (SPE)--> rappresenta la forma matematica che si vuole usare per il fit sul primo foto-elettrone
    gauss1 = par[0]*np.exp(-0.5*((x[0]-par[1])/par[2])**2)                         # ha tre parametri ( li inizializzo prima del fit globale sugli istogrammi su SPE )
    return float(gauss1)                                                           # il return float(gauss1) è necessario perchè.... same as above 


    


def double_gaus(x, par):                                                           # fit globale che si effettuerà sugli istogrammi di single-photo-electron (SPE)
    #return float(par[0]*np.exp(-0.5*((x[0]-par[1])/par[2])**2))
    #exp = par[0]*np.exp(-1*x[0]/par[1])
    gauss0 = par[0]*np.exp(-0.5*((x[0])/par[1])**2)                                # gaussiana di fit per il piedistallo ( o meglio per la coda destra della gaussiana del piedistallo che è visibile )
    gauss1 = par[2]*np.exp(-0.5*((x[0]-par[3])/par[4])**2)                         # gaussiana di fit per il single-photo-electron ( o meglio per la gaussiana del single-photo-electron  )
    gauss2 = par[5]*np.exp(-0.5*((x[0]-2*par[3])/par[4])**2)                       # gaussiana di fit per il secondo fotoelettrone ( o meglio per la gaussiana del secondo foto-elettrone )
    #return float(gauss0+gauss1+gauss2)
    return float(gauss0+gauss1)









def SPE():                                                                         # definisco la funzione SPE che prendendo i dati relativi al singolo fotoelettrone incidente ( condizione di bassa intensità, ovvero valor medio uguale ad 1, a cui comunque avrò dei conteggi a 0 e a 2 per via dei contegggi della poissoniana---> VALORE MEDIO UGUALE AD 1 NEL SENSO CHE AD OGNI IMPULSO DEL LED VIENE EMESSO IN MEDIA 1 FOTONE--> DISTRIBUZIONE DI POISSON PER EVENTI INDIPENDENTI, RARI E DISCRETI ) 
                                                                                   # esegue i vari fit sui dati in modo da inizializzare i parametri che poi saranno usati per fare il fit totale.
                                                                                   # Questo mi serve perchè in questa parte (lavoro inziale) siamo interessati a studiare il piedistallo di cui però visualizziamo solo la coda destra della gaussina corrispondente (gaussiana 0) e allora, per risolvere, facendo un fit che deve fittare bene gaussiana0+gaussiana1+gaussiana2 
                                                                                   # vogliamo di fatto, per ogni PMT, confrontare se il valore medio della distribuzione gaussiana del rumore dell'elettronica a PMT non alimentato sia uguale al valore medio della gaussiana 0 che si ricava facendo comparire il secondo picco guassiana e stimando la distanza fra i valori 
                                                                                   # medi delle gassiane 1 e 2 per poi inferire quello tra gaussiana 1 e 0 così da inferire il valor medio della guassiana 0 e confrontarlo, come detto, con i valori medi ottenuti ad elettronica spenta   

    fit_ped = ROOT.TF1("f", ped, 35, 80, 2)                                        # definisco la funzione di fit 'effettiva' per la gussiana 
                                                                                   # del piedistallo che 'incarna' la forma matematica che vogliamo 
                                                                                   # usare per il fit ( esponenziale ) e il range in cui viene fittata ( ovvero da 35 a 80 )
                                                                                   # e, quindi, la funzione che fittando la gaussiana del piedistallo inizializza i parametri
                                                                                   # che saranno poi usati come parametri del 'fit totale'.
                                                                                   # Sto creando un oggetto TF1 di ROOT, che rappresenta una funzione
                                                                                   # da usare per fare il fit su un istogramma

    fit_ped.SetParNames("A", "sigma0")                                             # sto dando i nomi ai parametri della funzione della funzione di fit_ped




    fit_prima = ROOT.TF1("s", primm, 70, 400, 3)                                   # definisco la funzione di fit 'effettiva' per la gaussiana 1 
                                                                                   # --> quella relativa al "single-photo-electron"
                                                                                   # che incarna la forma matematica che vogliamo usare per il fit 
                                                                                   # e il range in cui viene fittata ( ovvero fa 70 a 400 )
                                                                                   # sto creando un oggetto TF1 di ROOT, che rappresenta
                                                                                   # una funzione da usare per fare il fit su un istogramma

    fit_prima.SetParNames("B, mu, sigma")                                          # sto dando i nomi ai parametri della funzione della funzione di fit_prima


    #fit_function = ROOT.TF1("f", double_gaus, 35, 450, 6)
    #fit_function.SetParNames("A", "sigma0", "B", "mu", "sigma", "C")







    spe_data_long = Path("/swgo/multiPMT/calibration/batch_1/single_photoelectron/2025_03_17/run_spe_v_1450_thr_45_pol_90_whe_10_8_delay_7100_win_400")
    
    # sto definendo una variabile di nome spe_data_long che contiene il percorso alla cartella del file/ dei file che verranno letti 
    # esso è un oggetto della classe Path che è stata importata all'inizio del codice e su cui posso agire con i metodi che sono contenuti nella classe 
    # Path contenuta nella libreria pathlib
    # NOTA QUINDI CHE spe_data_long E' UN OGGETTO DI TIPO Path creato con la sintassi oggetto=NOME_CLASSE(ARGOMENTO)
    # la sintassi è questa perchè ho importato la classe all'interno dello script e quindi devo usare la sintassi indicata


    df_list = []                                                                     # Inizializza una lista vuota che conterrà i DataFrame letti dai file CSV


    for data_file in spe_data_long.rglob("*.csv"):                                   # .rglob è un metodo della classe Path (della libreria pathlib)
                                                                                     # Il metodo rglob("*.csv") restituisce un generatore di oggetti Path, uno per ciascun file CSV trovato ricorsivamente a partire dalla cartella collegata al percorso contenuto da spe_data_long.
                                                                                     # Un generatore di oggetti Path è un oggetto iterabile che produce i percorsi dei file uno alla volta, su richiesta. --> guarda note

                                                                                     # A QUESTO PROPOSITO NOTA CHE AVENDO IMPORTATO LA CLASSE Path DA pathlib ALL'INIZIO NON HO BISGNO DI SCRIVERE  Path.rglob(spe_data_long,"*.csv") ovvero classe.metodo(oggetto)
                                                                                     # MA INVECE POSSO SEMPLICEMENNTE SCRIVERE spe_data_long.rglob("*.csv") OVVERO oggetto.metodo
                                                                                     # Entrambe funzionano, ma in questo caso utilizziamo: istanza/oggetto.metodo() avendo importato la classe
                                                                                     # data_file ora è la variabile di iterazione del ciclo che ad ogni ciclo diventa il file, o meglio il percorso verso lo specifico file.
                                                                                     # Il metodo rglob("*.csv") restituisce un generatore di oggetti Path, uno per ciascun file CSV trovato ricorsivamente
                                                                                     # a partire dalla cartella spe_data_long.
                                                                                     # RICORDA: almeno devi importare il modulo o la libreria in cui è definita la classe, altrimenti Python non sa nemmeno che esiste.



        df = pd.read_csv(data_file)                                                  # legge il CSV e lo trasforma in un DataFrame ( tabella )
        df_list.append(df)                                                           # aggiunge/ appende alla lista df_list il dataframe appena creato dalla lettura
                                                                                     # del file letto a questa iterazione 

    df = pd.concat(df_list)                                                          # unisce tutti i DataFrame presenti nella lista di dataframe df_list in un unico DataFrame "finale", con nome df, che è l'unione di tutti i DataFrame (tabelle) che sono stati costruiti a partire dal singoli file di dati

    hist_list = []                                                                   # creo una lista vuota ( che sarà quella in cui alloco gli istogrammi )
    canvas = ROOT.TCanvas("a", "Multi-Channel Histogram", 1800, 1000)                # creo una canvas in cui visualizzerò gli istogrammi per i singoli PMT
    canvas.Divide(3, 3)

    guadagno_flt=[]                                                                   # definisco un array con i valori dei guadagni --> in verità questo array conterrà solo i valori di energia in pC e non direttamente i guadagni
    guadagno = []                                                                     # definisco un array con i valori dei guadagni 



    for channel in CHANNELS:                                                         # 'ciclo for' sul range CHANNELS (0-->6): questo vuol dire che la variabile Channel va da 0 a 6

        hist = ROOT.TH1D(f"hist_ch{channel}", f"Channel {channel}", 500, 0, 1000)    # viene creato l'istogramma, uno per ogni ciclo --> in questo senso si ha che ad ogni ciclo avrò
                                                                                     # un istogramma nuovo che come vedremo sarà riempito con i valori di energia relativi al PMT corrispondente 

        hist_list.append(hist)                                                       # l'istogramma viene aggiunto alla lista che abbiamo creato

        energy_values = df[df["Channel"] == channel]["Energy"].values                # all'interno della i-esima iterazione viene creato ora un array con i valori di energia contenuti nel DataFrame 'df' (tabella) corrispondenti 
                                                                                     # al PMT considerato alla iterazione i-esima.
                                                                                     # I valori inseriti nell'array alla i-esima iterazione del 'ciclo for' sono quindi i valori contenuti DataFrame 'df' (tabella), che contiene i valori di energia dei diversi PMT,
                                                                                     # che sono corrispondenti al PMT considerato in quella stessa iterazione ---> i diversi PMT ( numerati da 0 a 6 ) sono individuati dal valore assunto dalla variabile di iterazione 
                                                                                     # Channel che varia appunto da 0 a 6; valori di energia presi appunto dal DataFrame 'df'
                                                                                     # In questo senso ad ogni iterazione ( con Channel che assume valori da 0 a 6 ) trovo che i diversi valori di energia
                                                                                     # rivelati separatamente dai singoli PMT vanno a riempire il corrispondente istogramma ad ogni ciclo successivo


        for ener in energy_values:                                                  # 'ciclo for' "annidato" (perchè è dentro il ciclo for più esterno) sui valori di energia riferiti ad un singolo canale/PMT: sto di fatto ciclando sui valori di 
                                                                                    # energia che ho catturato con il primo ciclo for alla riga 202 ovvero estraendo i valori di energia dal dataframe 'df'

            hist.Fill(ener-pedestal[channel])                                       # riempio l'istogramma con i valori di energia dell'array "energy_values" che viene sovrascritto ad ongi ciclo e che 
                                                                                    # contiene per ogni ciclo i valori di energia che vengono da un singolo PMT n-esimo, per una iterazione i-esima del ciclo for iniziale

        print("Fit channel ",channel)

        
        canvas.cd(channel+1)
        hist.SetLineColor(ROOT.kBlue)
        hist.SetFillColorAlpha(ROOT.kBlue, 0.3)

        max_scale = int(hist.GetEntries()/10)

        ##################################
        fit_ped.SetParameters(max_scale*2, 40)                                     # setto i parametri della funzione di fit 'ped' che fitta la gaussiana del piedistallo nel range opportuno considerato 
        fit_ped.SetParLimits(0, 0, max_scale*5) #A ped                             # setto i limiti per il parametro 0 della funzione di fit 'ped' che fitta la gaussiana del piedistallo nel range opportuno considerato
        fit_ped.SetParLimits(1, 0, 60) #sigma                                      # setto i limiti per il parametro 1 della funzione di fit 'ped' che fitta la gaussiana del piedistallo nel range opportuno considerato
        hist.Fit(fit_ped, "RL")                                                    # fitto l'istogramma con la funzione di fit 'ped' che fitta la gaussiana del piedistallo nel range opportuno considerato
        #hist.Draw()
        #################################

        ##############################
        n_ped = fit_ped.GetParameter(0)                                             # estraggo il parametro 0 della funzione di fit 'ped' che fitta la gaussiana del piedistallo nel range opportuno considerato 
        sigma_ped = fit_ped.GetParameter(1)                                         # estraggo il parametro 1 della funzione di fit 'ped' che fitta la gaussiana del piedistallo nel range opportuno considerato 
        ####################################

        # print(n_ped, " ", sigma_ped)
        #print()
        #print("Double Gauss")



        fit_prima.SetParameters(max_scale*2, 300, 40)                               # setto i parametri della funzione di fit 'prima' che fitta la gaussiana del piedistallo nel range opportuno considerato 
    #   fit_prima.SetParLimits(0, 0, max_scale*5) #A ped
    #   fit_prima.SetParLimits(1, 0, 880) #mu
    #   fit_prima.SetParLimits(2, 0, 60) #sigma

        hist.Fit(fit_prima, "RL")                                                   # fitto l'istogramma con la funzione di fit 'ped' che fitta la gaussiana del piedistallo nel range opportuno considerato

        n_prima = fit_prima.GetParameter(0)                                           # estraggo il parametro 0 della gaussiana di singolo fotoelettrone
        mu_prima = fit_prima.GetParameter(1)                                          # estraggo il parametro 1 della gaussiana di singolo fotoelettrone
        sigma_prima = fit_prima.GetParameter(2)                                       # estraggo il parametro 1 della gaussiana di singolo fotoelettrone

        if mu_prima <= 200:                                                             # visto che come avevamo visto come per certi valori di mu serviva fittare la funzione totale in un certo range di parametri
                                                                                        # per far convergere il fit mentre per altri valori di mu occorreva fittare la funzione totale in un certo altro range di parametri
                                                                                        # poniamo un if per dividere questre possibilità una dalll'altra.

            fit_function = ROOT.TF1("f", double_gaus, 35, 300, 6)                       # Definiamo quindi di la funzione di fit totale 'effettiva' che 'incarna' la forma matematica che vogliamo usare per il fit
            fit_function.SetParNames("A", "sigma0", "B", "mu", "sigma", "C")            # e il RANGE in cui viene fittato sto creando un oggetto TF1 di ROOT, che rappresenta
                                                                                        # una funzione da usare per fare il fit su un istogramma
        
            opt = ROOT.Fit.DataOptions()
            rangeB = ROOT.Fit.DataRange()
            rangeB.SetRange(35, 300)
        else:                                                                          # altra possibilità che entra quando mu => 200
            fit_function = ROOT.TF1("f", double_gaus, 35, 450, 6)
            fit_function.SetParNames("A", "sigma0", "B", "mu", "sigma", "C")

        
        #fit_function.SetParameters(n_ped, sigma_ped, 2000, 200)
        fit_function.SetParameters(n_ped, sigma_ped, n_prima, mu_prima, sigma_prima, 100)   # sto facendo il fit 'totale' ponendo come parametri della funzione di fit totale PROPRIO QUELLI
                                                                                            # che avevo estratto dai fit singoli sulla gaussiana 0 e sulla gaussiana 1 separatamente 
        #fit_function.FixParameter(0, n_ped)
        #fit_function.FixParameter(1, sigma_ped)
        
    #    fit_function.SetParLimits(0, 0.8*n_ped, 1.2*n_ped) #A ped
    #    fit_function.SetParLimits(1, 0.8*sigma_ped, 1.2*sigma_ped) #l sigma
        
        fit_function.SetParLimits(2, 0.7*n_prima, 1.3*n_prima) #B prima gauss                # impongo dei limiti sui parametri che entrano nella funzione di fit 'totale'
        fit_function.SetParLimits(3, 0.7*mu_prima, 1.3*mu_prima) #mu 
        fit_function.SetParLimits(4, 0.7*sigma_prima, 1.3*sigma_prima) #sigma

    #    fit_function.SetParLimits(5, 0, max_scale/4) #C seconda gauss
        
        hist.Fit(fit_function, "R")                                                          # fitto l'istogramma con la funzione di fit 'totale'

        #hist.Draw("")
        #fit_ped.Draw("same")
        #fit_function.Draw("same")
        guadagno_flt.append((fit_function.GetParameter(3) * CALIBRAZIONE_ADC))               # calcolo il valore della carica di single-photo-electron in pC e lo inserisco nell'array guadagno_flt che pur chiamandosi 'guadagno'
                                                                                             # in verità contiene le cariche espresse in pC dei valori medi estratti delle gaussiane che fittano i single-photo-electron resitutiti dal fit 
                                                                                             # 'totale' una volta che sono settati come parametri di tale fit quelli restutiti dai fit fatti singolarmente

        valore = format((fit_function.GetParameter(3) * CALIBRAZIONE_ADC)/ CARICA_ELETTRONE, ".2e")   # calcolo il guadagno e lo inserisco nell'array contenente i guadagni per i diversi PMT che ora sono posti ad una 
                                                                                                      # tensione ben fissata ( 1450 V = tensione massima )
        guadagno.append(valore)

    
    print(guadagno)  
    canvas.Modified() 
    canvas.Update()
    print("Press any key to stop the programm...")
    input()

    return guadagno_flt






# IN QUESTA SECONDA PARTE CI PONIAMO AD UN INTENSITA' TALMENTE BASSA ( MINIMA INTENSITA') CHE IL PICCO GAUSSIANO DEL SECONDO FOTOELETTRONE E' TRASCURABILE E 
# QUINDI POSSO FITTARE SOLO UNA GAUSSIANA PER IL PICCO DELLO SPE INVECE CHE RICONSIDERARE NUOVAMENTE IL MODELLO CON IL 'fit totale' CHE ERA SERVITO PER CONFRONTARE 
# I VALORI MEDI DEL PIEDISTALLO PRESI A PMT NON ALIMENTATI CON QUELLI OTTENUTI DURANTE LA PRESA DATI.
# QUINDI IN QIESTE CONDIZIONI DI INTENSITA' DOPO AVER TRACCIATO GLI ISTOGRAMMI QUELLO CHE DEVO FARE E' CORRETTAMENTE FITTARE DELLE GAUSSIANE VISTO CHE ASSUMO CHE 
# NON VI SIANO ERRORI SISTEMATICI DURANTE LA PRESA DATI E CHE COMUNQUE IL PICCO DEL 2SPE SIA TRASCURABILE.

# In verita ciò che è scritto sopra è sbagliato visto che dopo aver carattierzzato il segnale di singolo fotoelettrone ciò che si fa è aumentare molto la intensità
# del fascio e porre i diversi PMT a diverse tensioni in modo da caratterizzare le curve di guadagno per i diversi PMT posti a diverse tensioni
# In questo senso avendo aumentato così tanto la intensità del fascio quello che accade è che la poissoniana per N grande tende alla distribuzione gaussiana e quindi vedrò delle gaussiane a valori grandi di energia dopo una acquisisizone dati

def GUADAGNO():

    info_gaudagno = {}             # crea il dizionario che come vedremo avrà la struttura { chiavi = voltaggi : elementi = (PMT n-esimo ; guadagno per quel PMT n-esimo posto a quel voltaggio, quello indicato dalla chaive)}  


    guadagno_data = Path("/swgo/multiPMT/calibration/batch_1/gain_curve/2025_03_20/run_gain_800_1450_thr_400_pol_75_whe_9_9_win_130")

    # sto definendo una variabile di nome guadagno_data che contiene il percorso alla cartella del file/ dei file che verranno letti 
    # esso è un oggetto della classe Path che è stata importata all'inizio del codice e su cui posso agire con i metodi che sono contenuti nella classe 
    # Path contenuta nella libreria pathlib
    # NOTA QUINDI CHE guadagno_data E' UN OGGETTO DI TIPO Path creato con la sintassi oggetto=NOME_CLASSE(ARGOMENTO)



    for data_file in guadagno_data.rglob("*.csv"):                                  # .rglob restituisce un iteratore che scorre tutti i file con estensione .csv ricorsivamente all'interno della cartella 'guadagno_data' e delle sue sottocartelle.
                                                                                    # Per ogni file trovato, il ciclo assegna quel file alla variabile data_file.
                                                                                    # Quindi, data_file sarà un oggetto Path (della libreria pathlib), che rappresenta il percorso completo del file CSV trovato in quel passaggio del ciclo.
        
        
                                                                                     # .rglob è un metodo della classe Path (classe della libreria pathlib)

                                                                                     # A QUESTO PROPOSITO NOTA CHE AVENDO IMPORTATO LA CLASSE Path DA pathlib ALL'INIZIO NON HO BISOGNO DI SCRIVERE  Path.rglob(guadagno_data,"*.csv") ovvero classe.metodo(oggetto)
                                                                                     # MA INVECE POSSO SEMPLICEMENNTE SCRIVERE spe_data_long.rglob("*.csv") OVVERO oggetto.metodo
                                                                                     # Entrambe funzionano, ma in questo caso utilizziamo: istanza.metodo() avendo importato la classe
                                                                                     # data_file ora è la variabile di iterazione del ciclo che ad ogni ciclo diventa il percorso al file con estensione .csv


        try:
            voltage = int(data_file.name.split("_")[-1].split(".")[0])               # estrae un numero (la tensione, in volt) dal nome di un file, supponendo che sia nel nome del file stesso questa informazione
        except:
            print("Impossibile recuperare il valore della tensione")
            return

        canvas = ROOT.TCanvas(f"{voltage}", f"Multi-Channel Histogram- V {voltage}", 1800, 1000) # creo una canvas dove disegnerò i singoli grafici di guadagno per ogni PMT singolo
        canvas.Divide(3, 3)

        hist_list = []                                                  # crea una lista di istogrammi
        print(voltage)                                                  # stampo il voltaggio che ho letto dalla stringa del nome del file 
        
        df = pd.read_csv(data_file)                                     # Legge il CSV e lo trasforma in un DataFrame ( tabella ) di nome df 
                                                                        # NOTA BENE CHE IL FILE .CSV CHE LEGGO E' IN GENERALE CONTIENE I DATI ACQUISITI PER TUTTI I PMT AD UN CERTO VOLTAGGIO BEN FISSATO: IN QUESTO CASO AVRO' PROPRIO CHE I DATI SONO I SEGNALI IN ENERGIA ( IN CANALI ADC ) CHE RIVELO COME OUTPUT DAI VARI PMT NEL TEMPO DELL'ACQUISIZIONE/ TEMPO DI PRESA DATI
                                                                        # QUINDI OGNI FILE RAPPRESENTA UNA CERTA ACQUISIZIONE AD UNA CERTA TENSIONE E CONTIENE TUTTI GLI "EVENTI" OVVERO TUTTI I SEGNALI IN USCITA CHE SONO STATI RIVEALATI IN USCITA DAI DIVERSI PMT. 
                                                                        # IN PARTICOLARE OGNI FOTONE PUO' GENERARE AL PIU' UN EVENTO NEL CASO IN CUI INTERAGISCA CON UN ELETTRONE DEL FOTOCATODO, SCALZANDOLO.  

        for channel in CHANNELS:                                        # 'ciclo for' sul range CHANNELS ( 0 --> 6 ): questo vuol dire che la variabile Channel va da 0 a 6

            hist = ROOT.TH1D(f"hist_ch{channel}_v_{voltage}", f"Channel {channel} ", 1000, 0, 16000)            # viene creato l'istogramma, uno per ogni ciclo --> in questo senso si ha che ad ogni ciclo avrò
                                                                                                                # un istogramma nuovo che come vedremo sarà riempito con i valori di energia relativi al PMT corrispondente ( TUTTO QUESTO AD UN CERTO VOLTAGGIO BEN DEFINITO --> DIPENDE SE L'INTERO CHE INDICA IL VOLTAGGIO E' STATO TROVATO O MENO )



            hist_list.append(hist)                                       # viene aggiunto l'istogramma alla lista di istogrammi prima creata --> è la lista di istogrammi per i vari PMT che si ha per la presa dati ad un certo voltaggio

            channel_data = df[df["Channel"] == channel]["Energy"]        # channel_data diventa un DataFrame ad una colonna con i valori di energia per il PMT i-esimo che corrisponde all'interazione i-esima del 'ciclo for'

            if not channel_data.empty:                                   # se channel data non è vuoto allora posso printare la tabella/ dataframe ad una colonna
                #print(channel_data)

                for ener in channel_data.values:                         # il metodo trasforma il dataframe ad una colonna in un array contenete i valori di energia per il PMT che stiamo considerando a quel ciclo i-esimo
                    hist.Fill(ener-pedestal[channel])                    # e adesso riempio l'istogramma creato a questa iterazione del ciclo for con quei valori di energia corrispondenti
                
                mean = hist.GetMean()                                    # prendo la media dei valori 
                sigma = hist.GetRMS()                                    # prendo la deviazione standard dei valori 
                amplitude = int(hist.GetEntries()/10)
                min_value = mean - 5*sigma                               # definisco il valore minimo del range del fit della gaussiana sugli istogrammi ( CONDIZIONE DI MINIMA INTENSITA' E PICCO DEL 2SPE TRASCURABILE )
                max_value = mean + 5*sigma                               # definisco il valore massimo del range del fit della gaussina sugli istogrammi ( CONDIZIONE DI MINIMA INTENSITA' E PICCO DEL 2SPE TRASCURABILE )
                hist.SetAxisRange(min_value, max_value, "X")

                fit_gauss = ROOT.TF1("f", "gaus", min_value, max_value, 3)      # fitto gli istogrammi con una guassina in un range fissato e dipendente dalla sigma 
                fit_gauss.SetParameters(amplitude, mean, sigma)                 # setto i parametri della gaussiana con cui vado a fittare i dati sperimentali 

                if voltage not in info_gaudagno:             # avendo creato un dizionario in cui vedo come chiavi le tensioni e poi vedo come elementi le tuple contenenti {il PMT i-esimo, il suo guadagno quella tensione specifica}

                    info_gaudagno[voltage] = []              # nella riga sopra si controlla che il valore di voltaggio non sia già presente come chiave nel dizionario e si esegue l'operazione di creazione della struttura { chiavi = voltaggi : elementi = (PMT n-esimo ; guadagno per quel PMT n-esimo posto a quel voltaggio, quello indicato dalla chaive)}
                                                             # Inoltre se la chiave è già presente e quindi se il valore di voltaggio è già presente allora non viene creata la struttura {chiave (voltaggio) ;(elementi....) }
                                                             # Alla fine di tutte le iterazioni sui diversi PMT quello che trovo è quindi una struttra del tipo  { chiavi = voltaggi : elementi = (PMT n-esimo ; guadagno per quel PMT n-esimo posto a quel voltaggio, quello indicato dalla chiave)}

                info_gaudagno[voltage].append((channel, fit_gauss.GetParameter(1)* CALIBRAZIONE_ADC))       # quindi si appende la tupla ( numero del canale, il valore del guadagno a una certa chiave = voltaggio) a quel voltaggio che fa come chiave e quindi poichè per un singolo file di dati ad un ben specifico voltaggio si fa tutto il processo di isolamento
                                                                                                            # dei dati per gli specifici PMT a quello specifico voltaggio vale che alla fine del ciclo sui canali per uno specifico voltaggio avrò che per quel voltaggio ci specifico il cui file ho analizzato vi saranno tutte le tuple (PMT n-esimo ; guadagno per quel
                                                                                                            # PMT n-esimo posto a quel voltaggio, quello indicato dalla chaive)  

                                                                                                            # alla fine di TUUTTOOO il processo E QUINDI DI TUTTI I CICLI SUI DIVERSI FILE RAPPRESENTANTI LE DIVERSE ACQUISIZIONI AI DIVERSI VOLTAGGI IN CUI, PER CIASCUNA, SI è CICLATO SUI DIVERSI CANALI PER COSTRUIRE I DIVERSI ISTOGRAMMI PER I DIVERSI PMT, avrò quindi che
                                                                                                            # ci sarà un dizionario con le chiavi=voltaggi e per ogni chiave le tuple contenenti { il PMT i-esimo, il suo guadagno quella tensione specifica }
                                                                                                            # si sta di fatto calcolando solo la carica e non il guadagno effettivamente)


            canvas.cd(channel+1)
            hist.SetLineColor(ROOT.kBlue)
            hist.SetFillColorAlpha(ROOT.kBlue, 0.3)
            hist.Fit(fit_gauss, "R")
            #hist.Draw("")

        canvas.Modified() 
        canvas.Update()
        print("Press any key to stop the programm...")
        input()
    
    return info_gaudagno




# QUESTA PARTE DEVE ESSERE INTERPRETATA A LIVELLO DI CODICE DOPO AVER ALZATO LA INTENSITA' LUMINOSA


def GAIN_CURVE():                                    # definisco la funzione per la curva curva di guadagno
    #ROOT.gROOT.SetBatch(True)
    carica_spe = SPE()                               # CARICA SPE() restituisce l'array per i vari PMT con le cariche uscenti come output dai singoli PMT dopo l'amplificazione di singolo fotoelettrone ( chiamata della funzione ) ALLA TENSIONE 1450 V
    carica_raw = GUADAGNO()                          # INTENSITA' AUMENTATA   # GUADAGNO() restituisce il dizionario con la struttura { chiavi = voltaggi : elementi = (PMT n-esimo ; guadagno per quel PMT n-esimo posto a quel voltaggio che però dato che non dividiamo per la carica dell'elettrone è ancora la carica in uscita dal PMT n-esimo IN CONDIZIONI DI ALTA LUMINOSITA')} ( chiamata della funzione )
    #ROOT.gROOT.SetBatch(False)


    npe = []                                         # definisco un array con i valori della carica rivelati una volta che abbiamo alzato la luminosità 
    tensioni_gaudagno = []                           # definisco un array con i valori delle tensioni 

    gain_curve = {}                                 # definisco un dizionario di nome gain_curve

    carica_1450 = carica_raw[1450]                  # definisco carica_1450 come l'array che contiene le tuple per V=1450V che sono del tipo
                                                    # ( channel = PMT, carica raccolta in output alla tensione 1450V per quel PMT [ sarebbe il guadagno se dividessimo per la carica dell'elettrone ma poichè non dividiamo troviamo solo la carica ] IN CONDIZIONI DI ALTA' LUMINOSITA' ) 

    #DOPO AVER ALZATO LA INTESNITà LUMINOSA FACCIO QUANTO E' SOTTO OVVERO VADO A CONSIDERARE PROPRIO L'ELEMENTO DEL DIZIONARIO CHE CONTIENE LE TUPLE CON (PMT n-esimo, CARICA RACCOLTA AVENDO AUMENTATO LA INTENSITA' E DIVIDO CIASCUNO DI ESSI PER IL VALORE DELLA CARICA DI SINGOLO FOTOELETTRONE TROVANDO COSI' IL NUMERO DI FOTOELETTRONI EMESSI DA CATODO E AMPLIFICATI DAL PMT, PER CIASCUN PMT, ALLE CONDIZIONI DI A CUI CI SIAMO POSTI--> INTENSITA' AUMENTATA MA STABILE)

    for data in carica_1450:                        # ciclo sugli elementi dell'array carica_1450 che ha come elemento le tuple ( channel=PMT, Guadagno relativo alla tensione 1450V per quel PMT )
        channel = data[0]                           # pongo come primo elemento di questo il PMT il canale
        charge = data[1]                            # pongo come secondo elemento dell'array la carica rivelata da questo PMT per il segnale generato da singolo fotoelettrone 
        npe_1450 = charge / (carica_spe[channel])   # ATTENZIONE CHE npe_1450 SONO IL NUMERO DI FOTOELETTRONI CHE VIENE EMESSO DAL CATODO IN SEGUITO ALL'INNALZAMENTO DELLA INTENSITA' LUMINOSA E calcola i valori di npe (number photo-electron) generati per ogni PMT considerando la carica raccolta dall n-esimo PMT IN CONDIZIONI DI ALTA' LUMINOSITA' e dividendola per la carica raccolta per il SPE al variare dei PMT, SEMPRE A QUEL VOLTAGGIO OVVERO A 1450 V, altrimenti non avrebbe senso tutto questo. Questo mi restituisce il numero di fotoni incidenti = fotoelettroni generati per ogni PMT ( se ammettiamo che ogni fotone generi un fotoelettrone). DA QUI CAPISCO CHE LA CARATTERIIZAZIONE DEL SEGNALE DI SPE E' FATTA ALLA TENSIONE MASSIMA 1450 V
        npe.append(npe_1450)                        # riempie l'array npe con i valori calcolati di npe_1450 ovvero con i numeri posti in successione dei fotoni incidenti sui singoli PMT FOTOELETTRONI


    for voltage in range(900, 1500, 50):            # a passi di 50 V vado da 900 V a 1500 V 

        carica = carica_raw[voltage]                # calcolo un array di tuple ciascuna indicante {il PMT n-esimo, la carica raccolta in output alla tensione tesione specifica nel ciclo corrente ( una valore da 900 v a 1500 v procedendo a passi di 50 V ) per quel PMT [ sarebbe il guadagno se dividessimo per la carica dell'elettrone ma poichè non dividiamo troviamo solo la carica ] }, per ogni singolo step di voltaggio 

        for data in carica:                         # ciclo sugli elementi dei canali e di nuovo

            channel = data[0]                       # ribattezzzo il primo dato delle tuple come channel --> sono i PMT
            charge = data[1]                        # ribattezzzo il secondo dato delle tuple come charge --> è la carica rivelata

            gaudagno = charge / (npe[channel] * CARICA_ELETTRONE)  # calcolo il guadagno, per un tensione fissata, di tutti i PMT

            if channel not in gain_curve:           
                gain_curve[channel] = []  

            gain_curve[channel].append((voltage, guadagno)) # ridefinsco la struttura a dizionario in cui ho i voltaggi e il guadagno per ogni canale/PMT 


    # FACCIO IN QUESTO MODO PERCHE' COSI' AVRO' CHE ANCHE SE DAI SINGOLI FOTOCATODI DEI PMT NON PARTE LO STESSO NUMERO DI FOTO ELETTRONI COMUNQUE HO CHE RIESCO TAMITE IL CICLO APPENA SOPRA A CALCOLARE I VALORI DEL GUADAGANO PER OGNI SINGOLO PMT AI DINVERSI VALORI DI TENSIONE 

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

        vertical_line = ROOT.TLine(intersezione, 1e6, intersezione, 1e7)
        vertical_lines.append(vertical_line)
        vertical_line.SetLineColor(ROOT.kGreen)
        vertical_line.SetLineWidth(2)
        vertical_line.Draw("same")
    
    canvas.Modified() 
    canvas.Update()
    print("Press any key to stop the programm...")
    input()

    





def PEDESTAL():                                                    # definisco la funzione 

    info_pedestal = {}

    pedestal_data = Path("dati_pedestal/")                         # creato il link da terminale associo a pedestal_data il path in cui si trova il file di dati dei piedistalli


    fit_ped = ROOT.TF1("f", primm, 875, 895, 3)                      # definisco la funzione di Fit che nota bene ha 3 parametri
    fit_ped.SetParNames("A", "mu", "sigma0")


   



    for data_file in pedestal_data.rglob("*.csv"):                 # ciclo sui file di dati presenti nella cartella che è individuata in maniera univoca dal path ( ho solo un file di dati in questo caso )
        df = pd.read_csv(data_file)                                # leggo il file e associo creo un dataframe (tabella) con i dati letti
        






    hist_list = []                                                               # creo l'array con gli istogrammi che riempiro nle ciclo for sui PMT
    canvas = ROOT.TCanvas("a", "Multi-Channel Histogram", 1800, 1000)            # creo la canvas
    canvas.Divide(3, 3)                                                          # divido la canvas

    



    canvas.Connect("Closed()", "TApplication", ROOT.gApplication, "Terminate()")


    for channel in CHANNELS:                                                         # ciclo sui diversi PMT

        hist = ROOT.TH1D(f"hist_ch{channel}", f"Channel {channel}", 30, 870, 900)    # creo l'istogramma relativo ai diversi PMT
        hist.GetYaxis().SetRangeUser(0, 2900)                                       

        hist_list.append(hist)                                                       # riempio la lista di istogrammi con gli istogrammi
        energy_values = df[df["Channel"] == channel]["Energy"].values


        for ener in energy_values:                                                   # ciclo for sui valori di energia presenti nell'array
            hist.Fill(ener)                                                          # riempio gli istogrammi con i valori di energia corrispondenti dei diversi isotgrammi
  
        print("Fit channel ",channel)

    
        canvas.cd(channel+1)
        hist.SetLineColor(ROOT.kBlue)
        hist.SetFillColorAlpha(ROOT.kBlue, 0.3)
    
        max_scale = int(hist.GetEntries()/10)
    
        fit_ped.SetParameters(2500, 885, 40)
        fit_ped.SetParLimits(0, 0, 3000) #A ped
        fit_ped.SetParLimits(1, 800, 900) #mu
        fit_ped.SetParLimits(2, 0, 100) #sigma
        #hist.Fit(fit_ped, "RL")
        #n_ped = fit_ped.GetParameter(0)
        #mu_ped = fit_ped.GetParameter(1)
        #sigma_ped = fit_ped.GetParameter(2)

        #print(n_ped, " ", mu_ped, " ",sigma_ped)
    
    
        hist.Fit(fit_ped, "R")                                                         # esco dal ciclo for che riempie l' istogramma creato all'iterazion e i-esima con i dati sperimentali relativi e  
        hist.Draw()                                                                     # disegno l'istogramma in quello specifico ciclo del ciclo for più esterno
        #fit_function.Draw("same")
    

    canvas.Modified() 
    canvas.Update()
    print("Press any key to stop the programm...")
    input()

    return info_pedestal








            


if __name__ == '__main__':

    #a = SPE()
    # print(a)


    SPE()
    GUADAGNO()
    PEDESTAL()
    GAIN_CURVE()
