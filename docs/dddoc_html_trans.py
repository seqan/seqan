def translate(text):

#spezielles
    text = text.replace("\\br", "<p>")   

#deutsche Umlaute und Sonderzeichen

    text = text.replace("\\\"a", "&auml;")
    text = text.replace("\\\"o", "&ouml;")
    text = text.replace("\\\"u", "&uuml;")
    text = text.replace("\\\"A", "&Auml;")
    text = text.replace("\\\"O", "&Ouml;")
    text = text.replace("\\\"U", "&Uuml;")
    text = text.replace("\\3", "&szlig;")
    text = text.replace("\\ss", "&szlig;")

#escapen spezieller Zeichen
    text = text.replace("\\colon", ":")
    text = text.replace("\\dot", ".")
    text = text.replace("\\at", "@")
    text = text.replace("\\pipe", "|")
    text = text.replace("\\dollar", "$")
    text = text.replace("\\quote", "\"")
    text = text.replace("\\backslash", "\\")

#Vergleiche

    text = text.replace("\\neq", "&#x2260;")
    text = text.replace("\\approx", "&#x2248;")
    text = text.replace("\\equiv", "&#x2261;")
    text = text.replace("\\leq", "&#x2264;")
    text = text.replace("\\geq", "&#x2265;")

#Binaere Operatoren

    text = text.replace("\\pm", "&plusmn;") #"<font face=\"symbol\">&#xb1;</font>")
    text = text.replace("\\times", "&times;") #"<font face=\"symbol\">&#xb4;</font>")
    text = text.replace("\\cdot", "&middot;") #"<font face=\"symbol\">&#xd7;</font>")
    text = text.replace("\\div", "&divide;") #"<font face=\"symbol\">&#xb8;</font>")
    text = text.replace("\\ast", "*")
    text = text.replace("\\circ", "&#x25cb;") #"<font face=\"symbol\">&#xb0;</font>")
    text = text.replace("\\otimes", "<font face=\"symbol\">&#xc4;</font>")
    text = text.replace("\\oplus", "<font face=\"symbol\">&#xc5;</font>")

#Logik

    text = text.replace("\\exists", "<font face=\"symbol\">&#x24;</font>")
    text = text.replace("\\forall", "<font face=\"symbol\">&#x22;</font>")
    text = text.replace("\\neg", "<font face=\"symbol\">&#xd8;</font>")
    text = text.replace("\\vee", "<font face=\"symbol\">&#xda;</font>")
    text = text.replace("\\wedge", "<font face=\"symbol\">&#xd9;</font>")

#Mengenlehre

    text = text.replace("\\in", "<font face=\"symbol\">&#xce;</font>")
    text = text.replace("\\ni", "<font face=\"symbol\">&#x27;</font>")
    text = text.replace("\\notin", "<font face=\"symbol\">&#xcf;</font>")
    text = text.replace("\\cup", "<font face=\"symbol\">&#xc8;</font>")
    text = text.replace("\\cap", "<font face=\"symbol\">&#xc7;</font>")
    text = text.replace("\\subset", "<font face=\"symbol\">&#xcc;</font>")
    text = text.replace("\\supset", "<font face=\"symbol\">&#xc9;</font>")
    text = text.replace("\\subseteq", "<font face=\"symbol\">&#xcd;</font>")
    text = text.replace("\\supseteq", "<font face=\"symbol\">&#xca;</font>")

#Pfeile

    text = text.replace("\\leftarrow", "&#x2190;")
    text = text.replace("\\rightarrow", "&#x2192;")
    text = text.replace("\\leftrightarrow", "&#x2194;")
    text = text.replace("\\Leftarrow", "<font face=\"symbol\">&#xdc;</font>")
    text = text.replace("\\Rightarrow", "<font face=\"symbol\">&#xde;</font>")
    text = text.replace("\\Leftrightarrow", "<font face=\"symbol\">&#xdb;</font>")

#Spezielle Zeichen

    text = text.replace("\\infty", "&#x221E;")
    text = text.replace("\\ldots", "...")
    text = text.replace("\\squared", "&#x00B2;")
    text = text.replace("\\cubic", "&#x00B3;")

#Griechische Buchstaben

    text = text.replace("\\Gamma", "&#x0393;")
    text = text.replace("\\Delta", "&#x0394;")
    text = text.replace("\\Theta", "&#x0398;")
    text = text.replace("\\Lambda", "&#x039b;")
    text = text.replace("\\Xi", "&#x039e;")
    text = text.replace("\\Pi", "&#x03a0;")
    text = text.replace("\\Sigma", "&#x03a3;")
    text = text.replace("\\Phi", "&#x03a6;")
    text = text.replace("\\Psi", "&#x03a8;")
    text = text.replace("\\Omega", "&#x03a9;")
    
    text = text.replace("\\alpha", "&#x03b1;")
    text = text.replace("\\beta", "&#x03b2;")
    text = text.replace("\\gamma", "&#x03b3;")
    text = text.replace("\\delta", "&#x03b4;")
    text = text.replace("\\epsilon", "&#x03b5;")
    text = text.replace("\\zeta", "&#x03b6;")
    text = text.replace("\\eta", "&#x03b7;")
    text = text.replace("\\theta", "&#x03b8;")
    text = text.replace("\\iota", "&#x03b9;")
    text = text.replace("\\kappa", "&#x03ba;")
    text = text.replace("\\lambda", "&#x03bb;")
    text = text.replace("\\mu", "&#x03bc;")
    text = text.replace("\\nu", "&#x03bd;")
    text = text.replace("\\xi", "&#x03be;")
    text = text.replace("\\omicron", "&#x03bf;")
    text = text.replace("\\pi", "&#x03c0;")
    text = text.replace("\\rho", "&#x03c1;")
    text = text.replace("\\varsigma", "&#x03c2;")
    text = text.replace("\\sigma", "&#x03c3;")
    text = text.replace("\\tau", "&#x03c4;")
    text = text.replace("\\upsilon", "&#x03c5;")
    text = text.replace("\\phi", "&#x03c6;")
    text = text.replace("\\chi", "&#x03c7;")
    text = text.replace("\\psi", "&#x03c8;")
    text = text.replace("\\omega", "&#x03c9;")

    return text
