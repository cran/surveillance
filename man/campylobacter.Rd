\name{campylobacter}
\alias{campylobacter}
\docType{data}
\title{Weekly Campylobacter Reports in Germany 2001 - (mid) 2003}
\description{
  Reported number of campylobacter cases for each "Kreis" in Germany
  over the period 1 January 2001 (IfSG) to June 2003. The data
  have been provided to us by the Robert Koch Institute (RKI), Berlin,
  Germany.
}
\usage{data(campylobacter)}
\format{
  A data frame with 130 weeks (rows) for all 439 Kreise (cols). Each
  entry is the number of reports in the Kreis during that week. See
  kreise.txt for information about the 439 kreise. Each "kreis" has a
  unique id -- in the data file this is the first line with
  headers. However these are not valid identifiers in R, why an "X" is
  prepended to the id number. 
  \describe{
\item{X1001}{SK Flensburg, which has a population of 84480}
\item{X1002}{SK Kiel, which has a population of 232242}
\item{X1003}{SK Lübeck, which has a population of 213496}
\item{X1004}{SK Neumünster, which has a population of 79646}
\item{X1051}{LK Dithmarschen, which has a population of 137447}
\item{X1053}{LK Herzogtum Lauenburg, which has a population of 181661}
\item{X1054}{LK Nordfriesland, which has a population of 165026}
\item{X1055}{LK Ostholstein, which has a population of 203386}
\item{X1056}{LK Pinneberg, which has a population of 293914}
\item{X1057}{LK Plön, which has a population of 133624}
\item{X1058}{LK Rendsburg-Eckernförde, which has a population of 271643}
\item{X1059}{LK Schleswig-Flensburg, which has a population of 198390}
\item{X1060}{LK Segeberg, which has a population of 252758}
\item{X1061}{LK Steinburg, which has a population of 136548}
\item{X1062}{LK Stormarn, which has a population of 219988}
\item{X2000}{SK Hamburg, which has a population of 1726363}
\item{X3101}{SK Braunschweig, which has a population of 245516}
\item{X3102}{SK Salzgitter, which has a population of 111696}
\item{X3103}{SK Wolfsburg, which has a population of 121887}
\item{X3151}{LK Gifhorn, which has a population of 172880}
\item{X3152}{LK Göttingen, which has a population of 265058}
\item{X3153}{LK Goslar, which has a population of 155610}
\item{X3154}{LK Helmstedt, which has a population of 99527}
\item{X3155}{LK Northeim, which has a population of 150100}
\item{X3156}{LK Osterode am Harz, which has a population of 84750}
\item{X3157}{LK Peine, which has a population of 132998}
\item{X3158}{LK Wolfenbüttel, which has a population of 126842}
\item{X3241}{Region Hannover, which has a population of 1122038}
\item{X3251}{LK Diepholz, which has a population of 212761}
\item{X3252}{LK Hameln-Pyrmont, which has a population of 162290}
\item{X3254}{LK Hildesheim, which has a population of 293091}
\item{X3255}{LK Holzminden, which has a population of 80866}
\item{X3256}{LK Nienburg (Weser), which has a population of 126310}
\item{X3257}{LK Schaumburg, which has a population of 166548}
\item{X3351}{LK Celle, which has a population of 181962}
\item{X3352}{LK Cuxhaven, which has a population of 206192}
\item{X3353}{LK Harburg, which has a population of 235110}
\item{X3354}{LK Lüchow-Dannenberg, which has a population of 51820}
\item{X3355}{LK Lüneburg, which has a population of 168734}
\item{X3356}{LK Osterholz, which has a population of 111473}
\item{X3357}{LK Rotenburg (Wümme), which has a population of 163058}
\item{X3358}{LK Soltau-Fallingbostel, which has a population of 141192}
\item{X3359}{LK Stade, which has a population of 192973}
\item{X3360}{LK Uelzen, which has a population of 97635}
\item{X3361}{LK Verden, which has a population of 133289}
\item{X3401}{SK Delmenhorst, which has a population of 76284}
\item{X3402}{SK Emden, which has a population of 51185}
\item{X3403}{SK Oldenburg, which has a population of 155908}
\item{X3404}{SK Osnabrück, which has a population of 164195}
\item{X3405}{SK Wilhelmshaven, which has a population of 84994}
\item{X3451}{LK Ammerland, which has a population of 111994}
\item{X3452}{LK Aurich, which has a population of 188363}
\item{X3453}{LK Cloppenburg, which has a population of 150973}
\item{X3454}{LK Emsland, which has a population of 304698}
\item{X3455}{LK Friesland, which has a population of 101402}
\item{X3456}{LK Grafschaft Bentheim, which has a population of 131086}
\item{X3457}{LK Leer, which has a population of 162765}
\item{X3458}{LK Oldenburg, which has a population of 122550}
\item{X3459}{LK Osnabrück, which has a population of 355871}
\item{X3460}{LK Vechta, which has a population of 128189}
\item{X3461}{LK Wesermarsch, which has a population of 94332}
\item{X3462}{LK Wittmund, which has a population of 57421}
\item{X4011}{SK Bremen, which has a population of 540950}
\item{X4012}{SK Bremerhaven, which has a population of 118701}
\item{X5111}{SK Düsseldorf, which has a population of 570765}
\item{X5112}{SK Duisburg, which has a population of 512030}
\item{X5113}{SK Essen, which has a population of 591889}
\item{X5114}{SK Krefeld, which has a population of 239559}
\item{X5116}{SK Mönchengladbach, which has a population of 262963}
\item{X5117}{SK Mülheim a.d.Ruhr, which has a population of 172332}
\item{X5119}{SK Oberhausen, which has a population of 221619}
\item{X5120}{SK Remscheid, which has a population of 118753}
\item{X5122}{SK Solingen, which has a population of 165032}
\item{X5124}{SK Wuppertal, which has a population of 364784}
\item{X5154}{LK Kleve, which has a population of 301861}
\item{X5158}{LK Mettmann, which has a population of 509023}
\item{X5162}{LK Neuss, which has a population of 445683}
\item{X5166}{LK Viersen, which has a population of 302795}
\item{X5170}{LK Wesel, which has a population of 476240}
\item{X5313}{SK Aachen, which has a population of 245778}
\item{X5314}{SK Bonn, which has a population of 306016}
\item{X5315}{SK Köln, which has a population of 967940}
\item{X5316}{SK Leverkusen, which has a population of 160829}
\item{X5354}{LK Aachen, which has a population of 307936}
\item{X5358}{LK Düren, which has a population of 271074}
\item{X5362}{LK Erftkreis, which has a population of 458430}
\item{X5366}{LK Euskirchen, which has a population of 190581}
\item{X5370}{LK Heinsberg, which has a population of 252275}
\item{X5374}{LK Oberbergischer Kreis, which has a population of 289518}
\item{X5378}{LK Rheinisch-Bergischer Kreis, which has a population of 276525}
\item{X5382}{LK Rhein-Sieg-Kreis, which has a population of 583051}
\item{X5512}{SK Bottrop, which has a population of 120780}
\item{X5513}{SK Gelsenkirchen, which has a population of 276740}
\item{X5515}{SK Münster, which has a population of 267197}
\item{X5554}{LK Borken, which has a population of 362834}
\item{X5558}{LK Coesfeld, which has a population of 216805}
\item{X5562}{LK Recklinghausen, which has a population of 656053}
\item{X5566}{LK Steinfurt, which has a population of 437545}
\item{X5570}{LK Warendorf, which has a population of 282286}
\item{X5711}{SK Bielefeld, which has a population of 323373}
\item{X5754}{LK Gütersloh, which has a population of 347954}
\item{X5758}{LK Herford, which has a population of 255393}
\item{X5762}{LK Höxter, which has a population of 155518}
\item{X5766}{LK Lippe, which has a population of 364990}
\item{X5770}{LK Minden-Lübbecke, which has a population of 323055}
\item{X5774}{LK Paderborn, which has a population of 292961}
\item{X5911}{SK Bochum, which has a population of 390087}
\item{X5913}{SK Dortmund, which has a population of 589240}
\item{X5914}{SK Hagen, which has a population of 202060}
\item{X5915}{SK Hamm, which has a population of 183805}
\item{X5916}{SK Herne, which has a population of 174018}
\item{X5954}{LK Ennepe-Ruhr-Kreis, which has a population of 349988}
\item{X5958}{LK Hochsauerlandkreis, which has a population of 281480}
\item{X5962}{LK Märkischer Kreis, which has a population of 457079}
\item{X5966}{LK Olpe, which has a population of 141705}
\item{X5970}{LK Siegen-Wittgenstein, which has a population of 295787}
\item{X5974}{LK Soest, which has a population of 307809}
\item{X5978}{LK Unna, which has a population of 430269}
\item{X6411}{SK Darmstadt, which has a population of 138457}
\item{X6412}{SK Frankfurt am Main, which has a population of 641076}
\item{X6413}{SK Offenbach, which has a population of 118429}
\item{X6414}{SK Wiesbaden, which has a population of 271076}
\item{X6431}{LK Bergstraße, which has a population of 264695}
\item{X6432}{LK Darmstadt-Dieburg, which has a population of 288463}
\item{X6433}{LK Groß-Gerau, which has a population of 250478}
\item{X6434}{LK Hochtaunuskreis, which has a population of 226683}
\item{X6435}{LK Main-Kinzig-Kreis, which has a population of 408144}
\item{X6436}{LK Main-Taunus-Kreis, which has a population of 221730}
\item{X6437}{LK Odenwaldkreis, which has a population of 100113}
\item{X6438}{LK Offenbach, which has a population of 336555}
\item{X6439}{LK Rheingau-Taunus-Kreis, which has a population of 184897}
\item{X6440}{LK Wetteraukreis, which has a population of 296056}
\item{X6531}{LK Gießen, which has a population of 254331}
\item{X6532}{LK Lahn-Dill-Kreis, which has a population of 263043}
\item{X6533}{LK Limburg-Weilburg, which has a population of 175679}
\item{X6534}{LK Marburg-Biedenkopf, which has a population of 253500}
\item{X6535}{LK Vogelsbergkreis, which has a population of 118475}
\item{X6611}{SK Kassel, which has a population of 194748}
\item{X6631}{LK Fulda, which has a population of 218326}
\item{X6632}{LK Hersfeld-Rotenburg, which has a population of 130010}
\item{X6633}{LK Kassel, which has a population of 245735}
\item{X6634}{LK Schwalm-Eder-Kreis, which has a population of 193404}
\item{X6635}{LK Waldeck-Frankenberg, which has a population of 170430}
\item{X6636}{LK Werra-Meißner-Kreis, which has a population of 113293}
\item{X7111}{SK Koblenz, which has a population of 107730}
\item{X7131}{LK Ahrweiler, which has a population of 130141}
\item{X7132}{LK Altenkirchen, which has a population of 137659}
\item{X7133}{LK Bad Kreuznach, which has a population of 157843}
\item{X7134}{LK Birkenfeld, which has a population of 90201}
\item{X7135}{LK Cochem-Zell, which has a population of 65860}
\item{X7137}{LK Mayen-Koblenz, which has a population of 212344}
\item{X7138}{LK Neuwied, which has a population of 185004}
\item{X7140}{LK Rhein-Hunsrück-Kreis, which has a population of 105591}
\item{X7141}{LK Rhein-Lahn-Kreis, which has a population of 129608}
\item{X7143}{LK Westerwaldkreis, which has a population of 202714}
\item{X7211}{SK Trier, which has a population of 100024}
\item{X7231}{LK Bernkastel-Wittlich, which has a population of 113998}
\item{X7232}{LK Bitburg-Prüm, which has a population of 96308}
\item{X7233}{LK Daun, which has a population of 64464}
\item{X7235}{LK Trier-Saarburg, which has a population of 138196}
\item{X7311}{SK Frankenthal, which has a population of 47827}
\item{X7312}{SK Kaiserslautern, which has a population of 99794}
\item{X7313}{SK Landau i.d.Pfalz, which has a population of 41089}
\item{X7314}{SK Ludwigshafen, which has a population of 162458}
\item{X7315}{SK Mainz, which has a population of 185293}
\item{X7316}{SK Neustadt a.d.Weinstraße, which has a population of 53822}
\item{X7317}{SK Pirmasens, which has a population of 44822}
\item{X7318}{SK Speyer, which has a population of 49956}
\item{X7319}{SK Worms, which has a population of 80625}
\item{X7320}{SK Zweibrücken, which has a population of 35799}
\item{X7331}{LK Alzey-Worms, which has a population of 125303}
\item{X7332}{LK Bad Dürkheim, which has a population of 134303}
\item{X7333}{LK Donnersbergkreis, which has a population of 78740}
\item{X7334}{LK Germersheim, which has a population of 124380}
\item{X7335}{LK Kaiserslautern, which has a population of 110207}
\item{X7336}{LK Kusel, which has a population of 78530}
\item{X7337}{LK Südliche Weinstraße, which has a population of 109694}
\item{X7338}{LK Ludwigshafen, which has a population of 147396}
\item{X7339}{LK Mainz-Bingen, which has a population of 196072}
\item{X7340}{LK Südwestpfalz, which has a population of 105271}
\item{X8111}{SK Stuttgart, which has a population of 587152}
\item{X8115}{LK Böblingen, which has a population of 367830}
\item{X8116}{LK Esslingen, which has a population of 505340}
\item{X8117}{LK Göppingen, which has a population of 257651}
\item{X8118}{LK Ludwigsburg, which has a population of 503229}
\item{X8119}{LK Rems-Murr-Kreis, which has a population of 412959}
\item{X8121}{SK Heilbronn, which has a population of 120163}
\item{X8125}{LK Heilbronn, which has a population of 324043}
\item{X8126}{LK Hohenlohekreis, which has a population of 108920}
\item{X8127}{LK Schwäbisch Hall, which has a population of 186967}
\item{X8128}{LK Main-Tauber-Kreis, which has a population of 137508}
\item{X8135}{LK Heidenheim, which has a population of 137055}
\item{X8136}{LK Ostalbkreis, which has a population of 315345}
\item{X8211}{SK Baden-Baden, which has a population of 53084}
\item{X8212}{SK Karlsruhe, which has a population of 279578}
\item{X8215}{LK Karlsruhe, which has a population of 422952}
\item{X8216}{LK Rastatt, which has a population of 224309}
\item{X8221}{SK Heidelberg, which has a population of 141509}
\item{X8222}{SK Mannheim, which has a population of 308385}
\item{X8225}{LK Neckar-Odenwald-Kreis, which has a population of 150091}
\item{X8226}{LK Rhein-Neckar-Kreis, which has a population of 527796}
\item{X8231}{SK Pforzheim, which has a population of 118002}
\item{X8235}{LK Calw, which has a population of 160115}
\item{X8236}{LK Enzkreis, which has a population of 193873}
\item{X8237}{LK Freudenstadt, which has a population of 121682}
\item{X8311}{SK Freiburg i.Breisgau, which has a population of 208294}
\item{X8315}{LK Breisgau-Hochschwarzwald, which has a population of 243043}
\item{X8316}{LK Emmendingen, which has a population of 152953}
\item{X8317}{LK Ortenaukreis, which has a population of 411491}
\item{X8325}{LK Rottweil, which has a population of 141761}
\item{X8326}{LK Schwarzwald-Baar-Kreis, which has a population of 211321}
\item{X8327}{LK Tuttlingen, which has a population of 133866}
\item{X8335}{LK Konstanz, which has a population of 268859}
\item{X8336}{LK Lörrach, which has a population of 219149}
\item{X8337}{LK Waldshut, which has a population of 166114}
\item{X8415}{LK Reutlingen, which has a population of 279177}
\item{X8416}{LK Tübingen, which has a population of 211058}
\item{X8417}{LK Zollernalbkreis, which has a population of 193196}
\item{X8421}{SK Ulm, which has a population of 118347}
\item{X8425}{LK Alb-Donau-Kreis, which has a population of 187000}
\item{X8426}{LK Biberach, which has a population of 184548}
\item{X8435}{LK Bodenseekreis, which has a population of 200948}
\item{X8436}{LK Ravensburg, which has a population of 270572}
\item{X8437}{LK Sigmaringen, which has a population of 133671}
\item{X9161}{SK Ingolstadt, which has a population of 117311}
\item{X9162}{SK München, which has a population of 1227958}
\item{X9163}{SK Rosenheim, which has a population of 59310}
\item{X9171}{LK Altötting, which has a population of 109023}
\item{X9172}{LK Berchtesgadener Land, which has a population of 100626}
\item{X9173}{LK Bad Tölz-Wolfratshausen, which has a population of 117416}
\item{X9174}{LK Dachau, which has a population of 131345}
\item{X9175}{LK Ebersberg, which has a population of 120416}
\item{X9176}{LK Eichstätt, which has a population of 121031}
\item{X9177}{LK Erding, which has a population of 117949}
\item{X9178}{LK Freising, which has a population of 155396}
\item{X9179}{LK Fürstenfeldbruck, which has a population of 195869}
\item{X9180}{LK Garmisch-Partenkirchen, which has a population of 87408}
\item{X9181}{LK Landsberg a.Lech, which has a population of 107609}
\item{X9182}{LK Miesbach, which has a population of 92720}
\item{X9183}{LK Mühldorf a.Inn, which has a population of 109733}
\item{X9184}{LK München, which has a population of 299883}
\item{X9185}{LK Neuburg-Schrobenhausen, which has a population of 90224}
\item{X9186}{LK Pfaffenhofen a.d.Ilm, which has a population of 113111}
\item{X9187}{LK Rosenheim, which has a population of 240106}
\item{X9188}{LK Starnberg, which has a population of 126588}
\item{X9189}{LK Traunstein, which has a population of 168827}
\item{X9190}{LK Weilheim-Schongau, which has a population of 128543}
\item{X9261}{SK Landshut, which has a population of 59513}
\item{X9262}{SK Passau, which has a population of 50669}
\item{X9263}{SK Straubing, which has a population of 44316}
\item{X9271}{LK Deggendorf, which has a population of 116558}
\item{X9272}{LK Freyung-Grafenau, which has a population of 82553}
\item{X9273}{LK Kelheim, which has a population of 110918}
\item{X9274}{LK Landshut, which has a population of 144504}
\item{X9275}{LK Passau, which has a population of 187355}
\item{X9276}{LK Regen, which has a population of 82721}
\item{X9277}{LK Rottal-Inn, which has a population of 119107}
\item{X9278}{LK Straubing-Bogen, which has a population of 96035}
\item{X9279}{LK Dingolfing-Landau, which has a population of 91218}
\item{X9361}{SK Amberg, which has a population of 44217}
\item{X9362}{SK Regensburg, which has a population of 127198}
\item{X9363}{SK Weiden i.d.OPf., which has a population of 43071}
\item{X9371}{LK Amberg-Sulzbach, which has a population of 109083}
\item{X9372}{LK Cham, which has a population of 131416}
\item{X9373}{LK Neumarkt i.d.OPf., which has a population of 127312}
\item{X9374}{LK Neustadt a.d.Waldnaab, which has a population of 101043}
\item{X9375}{LK Regensburg, which has a population of 178095}
\item{X9376}{LK Schwandorf, which has a population of 144370}
\item{X9377}{LK Tirschenreuth, which has a population of 79804}
\item{X9461}{SK Bamberg, which has a population of 69396}
\item{X9462}{SK Bayreuth, which has a population of 74519}
\item{X9463}{SK Coburg, which has a population of 42798}
\item{X9464}{SK Hof, which has a population of 50654}
\item{X9471}{LK Bamberg, which has a population of 143107}
\item{X9472}{LK Bayreuth, which has a population of 109249}
\item{X9473}{LK Coburg, which has a population of 92199}
\item{X9474}{LK Forchheim, which has a population of 112838}
\item{X9475}{LK Hof, which has a population of 108744}
\item{X9476}{LK Kronach, which has a population of 75458}
\item{X9477}{LK Kulmbach, which has a population of 78674}
\item{X9478}{LK Lichtenfels, which has a population of 70874}
\item{X9479}{LK Wunsiedel i.Fichtelgebirge, which has a population of 85278}
\item{X9561}{SK Ansbach, which has a population of 40412}
\item{X9562}{SK Erlangen, which has a population of 101912}
\item{X9563}{SK Fürth, which has a population of 111257}
\item{X9564}{SK Nürnberg, which has a population of 491307}
\item{X9565}{SK Schwabach, which has a population of 38518}
\item{X9571}{LK Ansbach, which has a population of 183668}
\item{X9572}{LK Erlangen-Höchstadt, which has a population of 129751}
\item{X9573}{LK Fürth, which has a population of 113530}
\item{X9574}{LK Nürnberger Land, which has a population of 168693}
\item{X9575}{LK Neustadt/Aisch-Bad Windsheim, which has a population of 99156}
\item{X9576}{LK Roth, which has a population of 124963}
\item{X9577}{LK Weißenburg-Gunzenhausen, which has a population of 95176}
\item{X9661}{SK Aschaffenburg, which has a population of 68361}
\item{X9662}{SK Schweinfurt, which has a population of 54568}
\item{X9663}{SK Würzburg, which has a population of 129915}
\item{X9671}{LK Aschaffenburg, which has a population of 174577}
\item{X9672}{LK Bad Kissingen, which has a population of 109498}
\item{X9673}{LK Rhön-Grabfeld, which has a population of 86717}
\item{X9674}{LK Haßberge, which has a population of 88412}
\item{X9675}{LK Kitzingen, which has a population of 89104}
\item{X9676}{LK Miltenberg, which has a population of 131369}
\item{X9677}{LK Main-Spessart, which has a population of 132206}
\item{X9678}{LK Schweinfurt, which has a population of 116874}
\item{X9679}{LK Würzburg, which has a population of 159311}
\item{X9761}{SK Augsburg, which has a population of 257836}
\item{X9762}{SK Kaufbeuren, which has a population of 42246}
\item{X9763}{SK Kempten, which has a population of 61505}
\item{X9764}{SK Memmingen, which has a population of 41101}
\item{X9771}{LK Aichach-Friedberg, which has a population of 124501}
\item{X9772}{LK Augsburg, which has a population of 238536}
\item{X9773}{LK Dillingen a.d.Donau, which has a population of 94435}
\item{X9774}{LK Günzburg, which has a population of 122114}
\item{X9775}{LK Neu-Ulm, which has a population of 161110}
\item{X9776}{LK Lindau, which has a population of 77850}
\item{X9777}{LK Ostallgäu, which has a population of 132389}
\item{X9778}{LK Unterallgäu, which has a population of 134920}
\item{X9779}{LK Donau-Ries, which has a population of 130335}
\item{X9780}{LK Oberallgäu, which has a population of 148315}
\item{X10041}{LK Stadtverband Saarbrücken, which has a population of 349948}
\item{X10042}{LK Merzig-Wadern, which has a population of 106222}
\item{X10043}{LK Neunkirchen, which has a population of 146678}
\item{X10044}{LK Saarlouis, which has a population of 211681}
\item{X10045}{LK Saar-Pfalz-Kreis, which has a population of 156732}
\item{X10046}{LK Sankt Wendel, which has a population of 95209}
\item{X11000}{SK Berlin, which has a population of 3388434}
\item{X12051}{SK Brandenburg a.d.Havel, which has a population of 76351}
\item{X12052}{SK Cottbus, which has a population of 105954}
\item{X12053}{SK Frankfurt (Oder), which has a population of 70308}
\item{X12054}{SK Potsdam, which has a population of 130435}
\item{X12060}{LK Barnim, which has a population of 171490}
\item{X12061}{LK Dahme-Spreewald, which has a population of 159568}
\item{X12062}{LK Elbe-Elster, which has a population of 129066}
\item{X12063}{LK Havelland, which has a population of 150340}
\item{X12064}{LK Märkisch-Oderland, which has a population of 189634}
\item{X12065}{LK Oberhavel, which has a population of 194022}
\item{X12066}{LK Oberspreewald-Lausitz, which has a population of 141959}
\item{X12067}{LK Oder-Spree, which has a population of 195670}
\item{X12068}{LK Ostprignitz-Ruppin, which has a population of 111878}
\item{X12069}{LK Potsdam-Mittelmark, which has a population of 212230}
\item{X12070}{LK Prignitz, which has a population of 94015}
\item{X12071}{LK Spree-Neiße, which has a population of 151100}
\item{X12072}{LK Teltow-Fläming, which has a population of 160414}
\item{X12073}{LK Uckermark, which has a population of 148606}
\item{X13001}{SK Greifswald, which has a population of 53533}
\item{X13002}{SK Neubrandenburg, which has a population of 71723}
\item{X13003}{SK Rostock, which has a population of 198964}
\item{X13004}{SK Schwerin, which has a population of 99978}
\item{X13005}{SK Stralsund, which has a population of 59970}
\item{X13006}{SK Wismar, which has a population of 46544}
\item{X13051}{LK Bad Doberan, which has a population of 118719}
\item{X13052}{LK Demmin, which has a population of 92935}
\item{X13053}{LK Güstrow, which has a population of 111205}
\item{X13054}{LK Ludwigslust, which has a population of 131292}
\item{X13055}{LK Mecklenburg-Strelitz, which has a population of 87115}
\item{X13056}{LK Müritz, which has a population of 69544}
\item{X13057}{LK Nordvorpommern, which has a population of 117722}
\item{X13058}{LK Nordwestmecklenburg, which has a population of 121153}
\item{X13059}{LK Ostvorpommern, which has a population of 113623}
\item{X13060}{LK Parchim, which has a population of 107996}
\item{X13061}{LK Rügen, which has a population of 74402}
\item{X13062}{LK Uecker-Randow, which has a population of 83459}
\item{X14161}{SK Chemnitz, which has a population of 255798}
\item{X14166}{SK Plauen, which has a population of 71155}
\item{X14167}{SK Zwickau, which has a population of 101726}
\item{X14171}{LK Annaberg, which has a population of 87200}
\item{X14173}{LK Chemnitzer Land, which has a population of 139814}
\item{X14177}{LK Freiberg, which has a population of 151560}
\item{X14178}{LK Vogtlandkreis, which has a population of 200388}
\item{X14181}{LK Mittlerer Erzgebirgskreis, which has a population of 93483}
\item{X14182}{LK Mittweida, which has a population of 136684}
\item{X14188}{LK Stollberg, which has a population of 93041}
\item{X14191}{LK Aue-Schwarzenberg, which has a population of 138331}
\item{X14193}{LK Zwickauer Land, which has a population of 133741}
\item{X14262}{SK Dresden, which has a population of 478631}
\item{X14263}{SK Görlitz, which has a population of 60264}
\item{X14264}{SK Hoyerswerda, which has a population of 47917}
\item{X14272}{LK Bautzen, which has a population of 155487}
\item{X14280}{LK Meißen, which has a population of 152044}
\item{X14284}{LK Niederschles. Oberlausitzkreis, which has a population of 103469}
\item{X14285}{LK Riesa-Großenhain, which has a population of 120539}
\item{X14286}{LK Löbau-Zittau, which has a population of 152304}
\item{X14287}{LK Sächsische Schweiz, which has a population of 145373}
\item{X14290}{LK Weißeritzkreis, which has a population of 124763}
\item{X14292}{LK Kamenz, which has a population of 154954}
\item{X14365}{SK Leipzig, which has a population of 493052}
\item{X14374}{LK Delitzsch, which has a population of 127162}
\item{X14375}{LK Döbeln, which has a population of 76681}
\item{X14379}{LK Leipziger Land, which has a population of 152687}
\item{X14383}{LK Muldentalkreis, which has a population of 135459}
\item{X14389}{LK Torgau-Oschatz, which has a population of 100485}
\item{X15101}{SK Dessau, which has a population of 80965}
\item{X15151}{LK Anhalt-Zerbst, which has a population of 76886}
\item{X15153}{LK Bernburg, which has a population of 69318}
\item{X15154}{LK Bitterfeld, which has a population of 106959}
\item{X15159}{LK Köthen, which has a population of 69791}
\item{X15171}{LK Wittenberg, which has a population of 129689}
\item{X15202}{SK Halle (Saale), which has a population of 243045}
\item{X15256}{LK Burgenlandkreis, which has a population of 141233}
\item{X15260}{LK Mansfelder Land, which has a population of 106523}
\item{X15261}{LK Merseburg-Querfurt, which has a population of 134105}
\item{X15265}{LK Saalkreis, which has a population of 81812}
\item{X15266}{LK Sangerhausen, which has a population of 67108}
\item{X15268}{LK Weißenfels, which has a population of 77712}
\item{X15303}{SK Magdeburg, which has a population of 229755}
\item{X15352}{LK Aschersleben-Staßfurt, which has a population of 102250}
\item{X15355}{LK Bördekreis, which has a population of 78926}
\item{X15357}{LK Halberstadt, which has a population of 78637}
\item{X15358}{LK Jerichower Land, which has a population of 99994}
\item{X15362}{LK Ohre-Kreis, which has a population of 117640}
\item{X15363}{LK Stendal, which has a population of 138553}
\item{X15364}{LK Quedlinburg, which has a population of 77783}
\item{X15367}{LK Schönebeck, which has a population of 76326}
\item{X15369}{LK Wernigerode, which has a population of 95188}
\item{X15370}{LK Altmarkkreis Salzwedel, which has a population of 100428}
\item{X16051}{SK Erfurt, which has a population of 200126}
\item{X16052}{SK Gera, which has a population of 109926}
\item{X16053}{SK Jena, which has a population of 101157}
\item{X16054}{SK Suhl, which has a population of 46765}
\item{X16055}{SK Weimar, which has a population of 63522}
\item{X16056}{SK Eisenach, which has a population of 44242}
\item{X16061}{LK Eichsfeld, which has a population of 113427}
\item{X16062}{LK Nordhausen, which has a population of 97685}
\item{X16063}{LK Wartburgkreis, which has a population of 143646}
\item{X16064}{LK Unstrut-Hainich-Kreis, which has a population of 118446}
\item{X16065}{LK Kyffhäuserkreis, which has a population of 92983}
\item{X16066}{LK Schmalkalden-Meiningen, which has a population of 142488}
\item{X16067}{LK Gotha, which has a population of 147418}
\item{X16068}{LK Sömmerda, which has a population of 80323}
\item{X16069}{LK Hildburghausen, which has a population of 73246}
\item{X16070}{LK Ilm-Kreis, which has a population of 121040}
\item{X16071}{LK Weimarer Land, which has a population of 90905}
\item{X16072}{LK Sonneberg, which has a population of 67175}
\item{X16073}{LK Saalfeld-Rudolstadt, which has a population of 131015}
\item{X16074}{LK Saale-Holzland-Kreis, which has a population of 93404}
\item{X16075}{LK Saale-Orla-Kreis, which has a population of 97568}
\item{X16076}{LK Greiz, which has a population of 122459}
\item{X16077}{LK Altenburger Land, which has a population of 112421}
  }
}
\source{Data have kindly been provided to us by the Robert
  Koch-Institut, Germany, 2003.}
\examples{
data(campylobacter)
#Show the number of cases for entire Germany
plot(apply(campylobacter,MARGIN=1,sum),type="l",ylab="reports",xlab="week no")
}
\keyword{}
