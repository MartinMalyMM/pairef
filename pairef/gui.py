# encoding: utf-8
# TO DO: unicode error
import sys
import os
import signal
import webbrowser
import platform
try:
    from PyQt4.QtGui import *
    from PyQt4.QtCore import *  #
except ImportError:
    sys.stderr.write(
        "Module PyQt4 is required for the graphical interface of PAIREF. "
        "Please install it or execute a PAIREF job from command-line.\n"
        "Tip: PyQt4 is installed in ccp4-python.\n"
        "Aborting.\n")
    sys.exit(1)


class MyWindow(QWidget):
    """PAIREF graphical interface (in PyQt4)"""
    def __init__(self):
        super(MyWindow, self).__init__()
        # Allow to term the app pressing Ctrl+C in console
        signal.signal(signal.SIGINT, signal.SIG_DFL)
        # locale.setlocale(locale.LC_ALL, "en_US.UTF-8")  useless
        # QLocale.setDefault(QLocale(QLocale.English, QLocale.UnitedStates))
        self.myInit()

    def myInit(self):
        """Specifications of widgets, their layout and connection"""
        self.setWindowTitle("PAIREF")
        logofilename = \
            str(os.path.dirname(__file__)) + os.path.sep + \
            'static' + os.path.sep + 'pairef_logo_64.png'
        self.setStyleSheet("""
QPushButton#required {
    font-weight: bold;
}

QLineEdit#required {
    background-color: lightyellow;
}

#indent {
    margin-left: 2em;
}

#italic {
    font-style: italic;
}
""")
        intValidatorPositive = QIntValidator()  # only real numbers
        intValidatorPositive.setRange(1, 999999)
        doubleValidator = QDoubleValidator()  # only real numbers
        # doubleValidatorPositive = QDoubleValidator()  # only real positive numbers
        # doubleValidatorPositive.setRange(0.000001, 999999)  # Not working well
        doubleValidatorPositive = doubleValidator  # Accepts also negative :-(
        alphanumValidator = QRegExpValidator(QRegExp("[A-Za-z0-9_-]+"))
        
        if os.path.isfile(logofilename):
            self.logoLabel = QLabel(self)
            logo = QPixmap(logofilename)
            self.logoLabel.setPixmap(logo)
            self.logoLabel.setAlignment(Qt.AlignCenter)
        self.header = QLabel("Automatic PAIRed REFinement protocol", 
                             self, objectName="italic")
        self.header.setAlignment(Qt.AlignCenter)

        resolutionRadio = QButtonGroup(self)
        self.step_numberRadio = QRadioButton("Define resolution step")
        resolutionRadio.addButton(self.step_numberRadio)
        self.step_numberRadio.setChecked(True)
        self.step_numberRadio.toggled.connect(self.resToggled)
        self.stepLabel = QLabel("Add high-resolution shells with a width of",
                                self, objectName="indent")
        self.stepLabel.setToolTip("in angstrom")
        self.step = QLineEdit("0.05")
        self.step.setValidator(doubleValidatorPositive)
        self.step.setToolTip("in angstrom")
        self.numberLabel = QLabel("Number of high-resolution shells",
                                  self, objectName="indent")
        self.number = QSpinBox()
        self.number.setMinimum(1)
        self.number.setValue(1)

        explicitToolTip = \
            "values in angstrom must be divided using commas without any " \
            "spaces and written in decreasing order"
        self.explicitRadio = QRadioButton(
            "Explicit definition of high-resolution shells")
        resolutionRadio.addButton(self.explicitRadio)
        self.explicitRadio.toggled.connect(self.resToggled)
        self.explicit = QLineEdit()
        self.explicit.setToolTip(explicitToolTip)
        self.explicit.setEnabled(False)
        self.explicitLabel = QLabel("e.g. 2.1,2.0,1.9",
                                    self, objectName="italic")
        self.explicitLabel.setToolTip("in angstrom")
        self.explicitLabel.setEnabled(False)

        self.res_init = QLineEdit()
        self.res_init.setToolTip("in angstrom")
        self.res_init.setValidator(doubleValidatorPositive)
        self.res_initLabel = QLabel("Starting resolution", self)
        self.res_initLabel.setToolTip("in angstrom")

        self.refinementLabel = QLabel("Refinement software", self)
        refinementRadio = QButtonGroup(self)
        self.refmacRadio = QRadioButton("REFMAC5", objectName="italic")
        refinementRadio.addButton(self.refmacRadio)
        self.refmacRadio.setChecked(True)
        self.refmacRadio.toggled.connect(self.refinementToggled)
        self.phenixRadio = QRadioButton("phenix.refine", objectName="italic")
        refinementRadio.addButton(self.phenixRadio)
        self.phenixRadio.setChecked(False)
        self.phenixRadio.toggled.connect(self.refinementToggled)

        self.projectName = QLineEdit("project")
        self.projectName.setValidator(alphanumValidator)
        self.projectLabel = QLabel("Project name", self)

        self.workdirName = QLineEdit(os.getcwd())
        self.workdirName.textChanged.connect(self.setworkdir)
        self.workdirName.setStyleSheet("background-color: lightgreen;")
        self.workdirButton = QPushButton("Working directory")
        self.workdirButton.clicked.connect(self.getworkdir)
        
        files = {}
        files["xyzin"] = {"label": "Structure model",
                          "filetype": "*.pdb *.mmcif *.cif"}
        files["hklin"] = {"label": "Diffraction data (merged)",
                          "filetype": "*.mtz"}
        files["hklin_unmerged"] = {"label": "Diffraction data (unmerged)",
                                   "filetype": "*.*"}
        files["libin"] = {"label": "Geometric restraints",
                          "filetype": "*.cif *.geo"}
        files["comin"] = {"label": "Command file for REFMAC5",
                          "filetype": "*.*"}
        files["def"] = {
            "label": "Keyword file for phenix.refine",
            "filetype": "*.def"}
        files["tlsin"] = {"label": "Input fixed TLS parameters",
                          "filetype": "*.tls *.tlsin *.tlsout"}
        
        self.xyzinFileName = QLineEdit(objectName="required")
        self.xyzinFileName.textChanged.connect(self.checkfile)
        self.xyzinButton = QPushButton(
            files["xyzin"]["label"], objectName="required")
        self.xyzinButton.clicked.connect(
            lambda: self.getfile(self.xyzinFileName, attrs=files["xyzin"]))

        self.hklinFileName = QLineEdit(objectName="required")
        self.hklinFileName.textChanged.connect(self.checkfile)
        self.hklinButton = QPushButton(
            files["hklin"]["label"], objectName="required")
        self.hklinButton.clicked.connect(
            lambda: self.getfile(self.hklinFileName, attrs=files["hklin"]))

        self.hklin_unmergedFileName = QLineEdit()
        self.hklin_unmergedFileName.textChanged.connect(self.checkfile)
        self.hklin_unmergedButton = QPushButton(
            files["hklin_unmerged"]["label"])
        self.hklin_unmergedButton.clicked.connect(
            lambda: self.getfile(self.hklin_unmergedFileName, attrs=files["hklin_unmerged"]))

        self.libinFileName = QLineEdit()
        self.libinFileName.textChanged.connect(self.checkfile)
        self.libinButton = QPushButton(
            files["libin"]["label"])
        self.libinButton.clicked.connect(
            lambda: self.getfile(self.libinFileName, attrs=files["libin"]))

        self.cominFileName = QLineEdit()
        self.cominFileName.textChanged.connect(self.checkfile)
        self.cominButton = QPushButton(
            files["comin"]["label"])
        self.cominButton.clicked.connect(
            lambda: self.getfile(self.cominFileName, attrs=files["comin"]))

        self.defFileName = QLineEdit()
        self.defFileName.textChanged.connect(self.checkfile)
        self.defButton = QPushButton(
            files["def"]["label"])
        self.defButton.clicked.connect(
            lambda: self.getfile(self.defFileName, attrs=files["def"]))
        self.defFileName.hide()
        self.defButton.hide()

        self.ncyc = QSpinBox()
        self.ncyc.setMinimum(1)
        self.ncyc.setValue(10)
        self.ncycLabel = QLabel("Number of refinement cycles",
                                self)
        self.ncyc.setToolTip("Number of refinement cycles that will "
                             "be performed in every resolution step")
        self.ncycLabel.setToolTip("Number of refinement cycles that will "
                                  "be performed in every resolution step")

        self.weightCheckBox = QCheckBox("Use automatic weighting")
        self.weightCheckBox.setChecked(True)
        self.weightCheckBox.stateChanged.connect(self.weight_show_hide)

        self.weightLabel = QLabel("Use weighting term")
        self.weightLabel.hide()
        self.weight = QLineEdit(objectName="required")
        self.weight.setValidator(doubleValidatorPositive)
        self.weight.textChanged.connect(self.checkweight)
        self.weight.hide()

        # TLS
        self.tlsCheckBox = QCheckBox("Input fixed TLS parameters")
        self.tlsCheckBox.stateChanged.connect(self.tls_show_hide)

        self.tlsinFileName = QLineEdit()
        self.tlsinFileName.textChanged.connect(self.checkfile)
        self.tlsinButton = QPushButton(
            files["tlsin"]["label"], objectName="indent")
        self.tlsinButton.clicked.connect(
            lambda: self.getfile(self.tlsinFileName, attrs=files["tlsin"]))

        self.tls_ncyc = QSpinBox()
        self.tls_ncyc.setMinimum(1)
        self.tls_ncyc.setValue(10)
        self.tls_ncycLabel = QLabel("Number of cycles of TLS refinement", self,
                                    objectName="indent")

        self.free = QLineEdit("0")
        self.free.setValidator(intValidatorPositive)
        self.freeLabel = QLabel("Exclude reflections with free-flag", self)

        self.tlsin_keep = QCheckBox("Keep using the same TLS input file "
                                    "in all the refinement runs",
                                    objectName="indent")

        tls_widgets = [self.tlsinFileName, self.tlsinButton, self.tls_ncyc,
                       self.tls_ncycLabel, self.tlsin_keep]
        for tls_widget in tls_widgets:
            tls_widget.hide()

        # Complete cross-validation, prerefinement, input model modification
        self.completeCheckBox = QCheckBox(
            "Perform complete cross-validation accross all the free-flag sets")
        self.completeCheckBox.stateChanged.connect(
            self.complete_show_hide)

        self.prerefinementCheckBox = QCheckBox(
            "Refine the input structure model at the starting "
            "resolution")
        self.prerefinementCheckBox.stateChanged.connect(
            self.prerefinement_show_hide)

        prerefinement_ncycToolTip = \
            "Number of refinement cycles to be performed as pre-" + \
            "refinement of the input structure model before paired " + \
            "refinement (the initial high resolution limit is " + \
            "used). Pre-refinement is performed by default in case " + \
            "of the complete cross-validation protocol."
        self.prerefinement_ncyc = QSpinBox()
        self.prerefinement_ncyc.setMinimum(1)
        self.prerefinement_ncyc.setValue(20)
        self.prerefinement_ncycLabel = QLabel(
            "Number of pre-refinement cycles", self, objectName="indent")
        self.prerefinement_ncyc.setToolTip(prerefinement_ncycToolTip)
        self.prerefinement_ncycLabel.setToolTip(prerefinement_ncycToolTip)

        prerefinement_shakeToolTip = \
            "Randomize coordinates of the input structure model " + \
            "with the given mean error value. This is done by " + \
            "default in the case of the complete cross-validation protocol."
        self.prerefinement_shakeCheckBox = QCheckBox(
            "Randomize coordinates with mean error", objectName="indent")
        self.prerefinement_shakeCheckBox.setToolTip(prerefinement_shakeToolTip)
        self.prerefinement_shake = QLineEdit("0.25")
        self.prerefinement_shake.setValidator(doubleValidatorPositive)
        self.prerefinement_shake.setToolTip(prerefinement_shakeToolTip)

        self.prerefinement_resetCheckBox = QCheckBox(
            "Reset B-factors to mean value", objectName="indent")
        self.prerefinement_resetCheckBox.stateChanged.connect(
            self.prerefinement_setUncheck)
        self.prerefinement_resetCheckBox.setToolTip(
            "Reset atomic B-factors of the input structure model to "
            "the mean value. This is done by default in the case of "
            "the completecross-validation protocol.")

        prerefinement_setToolTip = \
            "Set atomic B-factors of the input structure model to "
        "the given value (angstrom)"
        self.prerefinement_setCheckBox = QCheckBox(
            "Set atomic B-factors", objectName="indent")
        self.prerefinement_setCheckBox.setToolTip(prerefinement_setToolTip)
        self.prerefinement_setCheckBox.stateChanged.connect(
            self.prerefinement_resetUncheck)
        self.prerefinement_set = QLineEdit("30.0")
        self.prerefinement_set.setValidator(doubleValidatorPositive)
        self.prerefinement_set.setToolTip(prerefinement_setToolTip)

        prerefinement_addToolTip = \
            "Add the given value to B-factors of the input " + \
            "structure model (in angstrom)"
        self.prerefinement_addCheckBox = QCheckBox(
            "Add value to B-factors", objectName="indent")
        self.prerefinement_addCheckBox.setToolTip(prerefinement_addToolTip)
        self.prerefinement_add = QLineEdit("0")
        self.prerefinement_add.setValidator(doubleValidator)
        self.prerefinement_add.setToolTip(prerefinement_addToolTip)

        prerefinement_widgets = [
            self.prerefinement_ncyc, self.prerefinement_ncycLabel, 
            self.prerefinement_addCheckBox, self.prerefinement_add,
            self.prerefinement_setCheckBox, self.prerefinement_set,
            self.prerefinement_resetCheckBox,
            self.prerefinement_shakeCheckBox, self.prerefinement_shake]
        for prerefinement_widget in prerefinement_widgets:
            prerefinement_widget.hide()


        self.runButton = QPushButton("RUN")
        self.runButton.clicked.connect(self.run)
        self.stopButton = QPushButton("STOP")
        self.stopButton.clicked.connect(self.stop)
        self.stopButton.hide()

        self.browserLabel = QLabel("Job is running, open HTML log file ")
        self.browserLabel.hide()
        self.browserButton = QPushButton("Open in browser...")
        self.browserButton.hide()
        self.browserButton.clicked.connect(self.openBrowser)
        # myButton.move(150, 80)
        
        # self.setGeometry(96, 56, 406, 200)
        self.resize(800, 240)
        box_head = QVBoxLayout()
        grid_res = QGridLayout()
        grid_files = QGridLayout()
        grid_more = QGridLayout()
        grid_more2 = QGridLayout()
        grid_foot = QGridLayout()
        # grid.setSpacing(4)
        if os.path.isfile(logofilename):
            box_head.addWidget(self.logoLabel)
        box_head.addWidget(self.header)
        box_head.addStretch()
        grid_res.addWidget(self.step_numberRadio, 150, 0)
        grid_res.addWidget(self.stepLabel, 152, 0)
        grid_res.addWidget(self.step, 152, 2)
        grid_res.addWidget(self.numberLabel, 152, 4)
        grid_res.addWidget(self.number, 152, 6)
        grid_res.addWidget(self.explicitRadio, 156, 0)
        grid_res.addWidget(self.explicit, 156, 2)
        grid_res.addWidget(self.explicitLabel, 156, 4)
        grid_res.addWidget(self.res_initLabel, 180, 0)
        grid_res.addWidget(self.res_init, 180, 2)
        grid_res.addWidget(self.refinementLabel, 190, 0)
        grid_res.addWidget(self.refmacRadio, 190, 2)
        grid_res.addWidget(self.phenixRadio, 190, 4)
        grid_files.addWidget(self.projectLabel, 200, 0)
        grid_files.addWidget(self.projectName, 200, 2)
        grid_files.addWidget(self.workdirButton, 250, 0)
        grid_files.addWidget(self.workdirName, 250, 2)
        grid_files.addWidget(self.xyzinButton, 300, 0)
        grid_files.addWidget(self.xyzinFileName, 300, 2)
        grid_files.addWidget(self.hklinButton, 400, 0)
        grid_files.addWidget(self.hklinFileName, 400, 2)
        grid_files.addWidget(self.hklin_unmergedButton, 500, 0)
        grid_files.addWidget(self.hklin_unmergedFileName, 500, 2)
        grid_files.addWidget(self.libinButton, 600, 0)
        grid_files.addWidget(self.libinFileName, 600, 2)
        grid_files.addWidget(self.cominButton, 700, 0)
        grid_files.addWidget(self.cominFileName, 700, 2)
        grid_files.addWidget(self.defButton, 700, 0)
        grid_files.addWidget(self.defFileName, 700, 2)

        grid_more.addWidget(self.ncycLabel, 730, 0)
        grid_more.addWidget(self.ncyc, 730, 2)
        grid_more.addWidget(self.freeLabel, 730, 4)
        grid_more.addWidget(self.free, 730, 6)

        grid_more2.addWidget(self.weightCheckBox, 736, 0)
        grid_more2.addWidget(self.weightLabel, 740, 0)
        grid_more2.addWidget(self.weight, 740, 2)

        grid_more2.addWidget(self.tlsCheckBox, 750, 0)
        grid_more2.addWidget(self.tlsinButton, 760, 0)
        grid_more2.addWidget(self.tlsinFileName, 760, 2)
        grid_more2.addWidget(self.tls_ncycLabel, 761, 0)
        grid_more2.addWidget(self.tls_ncyc, 761, 2)
        grid_more2.addWidget(self.tlsin_keep, 762, 0, 1, 3)

        grid_more2.addWidget(self.completeCheckBox, 950, 0, 1, 3)
        grid_more2.addWidget(self.prerefinementCheckBox, 955, 0, 1, 3)
        grid_more2.addWidget(self.prerefinement_ncycLabel, 960, 0)
        grid_more2.addWidget(self.prerefinement_ncyc, 960, 2)
        grid_more2.addWidget(self.prerefinement_shakeCheckBox, 962, 0)
        grid_more2.addWidget(self.prerefinement_shake, 962, 2)
        grid_more2.addWidget(self.prerefinement_resetCheckBox, 964, 0)
        grid_more2.addWidget(self.prerefinement_setCheckBox, 965, 0)
        grid_more2.addWidget(self.prerefinement_set, 965, 2)
        grid_more2.addWidget(self.prerefinement_addCheckBox, 968, 0)
        grid_more2.addWidget(self.prerefinement_add, 968, 2)

        grid_foot.addWidget(self.runButton, 1000, 0)
        grid_foot.addWidget(self.stopButton, 1000, 1)
        grid_foot.addWidget(self.browserLabel, 1100, 0)
        grid_foot.addWidget(self.browserButton, 1100, 1)

        layout = box_head
        layout.addLayout(grid_res)
        layout.addLayout(grid_files)
        layout.addLayout(grid_more)
        layout.addLayout(grid_more2)
        layout.addLayout(grid_foot)
        self.setLayout(layout)
        self.show()

    def getworkdir(self):
        workdirname = QFileDialog.getExistingDirectory(self, 'Select directory')
        self.workdirName.setText(workdirname)

    def setworkdir(self):
        try:
            if os.path.isdir(self.workdirName.text()):
                os.chdir(self.workdirName.text())
                self.workdirName.setStyleSheet("background-color: lightgreen;")
            else:
                self.workdirName.setStyleSheet("background-color: red;")
        except UnicodeDecodeError:
            self.workdirName.setStyleSheet("background-color: red;")
            
    def getfile(self, QLine, attrs):
        action = "Select " + attrs["label"].lower()
        filetype = attrs["label"] + " (" + attrs["filetype"] + ")" + \
            " ;; All files (*.*)"
        if "win" in platform.system().lower() or \
                "mac" in platform.system().lower():
            # file type filter does not work well
            filename = QFileDialog.getOpenFileName(
                self, action, "")
        else:
            filename = QFileDialog.getOpenFileName(
                self, action, "", filetype)
        QLine.setText(filename)

    def checkfile(self):
        try:
            if os.path.isfile(self.sender().text()):
                self.sender().setStyleSheet("background-color: lightgreen;")
            else:
                self.sender().setStyleSheet("background-color: red;")
        except UnicodeDecodeError:
            self.sender().setStyleSheet("background-color: red;")

    def checkweight(self):
        weight = str(self.weight.text())
        try:
            if float(weight) > 0:
                self.sender().setStyleSheet("background-color: lightgreen;")
            else:
                self.sender().setStyleSheet("background-color: red;")
        except:
            self.sender().setStyleSheet("background-color: red;")

    def refinementToggled(self):
        widgets_refmac = [
            self.cominButton, self.cominFileName, self.weightCheckBox,
            self.tlsCheckBox]
        widgets_refmac_weight = [self.weightLabel, self.weight]
        widgets_refmac_tls = [
            self.tlsinButton, self.tlsinFileName, self.tls_ncycLabel,
            self.tlsin_keep, self.tls_ncyc]
        widgets_phenix = [self.defButton, self.defFileName]
        if self.refmacRadio.isChecked():
            self.ncyc.setValue(10)
            if not self.weightCheckBox.isChecked():
                widgets_refmac += widgets_refmac_weight
            if self.tlsCheckBox.isChecked():
                widgets_refmac += widgets_refmac_tls
            for widget in widgets_refmac:
                widget.show()
            for widget in widgets_phenix:
                widget.hide()
        if self.phenixRadio.isChecked():
            self.ncyc.setValue(3)
            widgets_refmac += widgets_refmac_weight + widgets_refmac_tls
            for widget in widgets_refmac:
                widget.hide()
            for widget in widgets_phenix:
                widget.show()
            

    def resToggled(self):
        if self.step_numberRadio.isChecked():
            self.step.setEnabled(True)
            self.number.setEnabled(True)
            self.explicit.setEnabled(False)
            self.explicitLabel.setEnabled(False)
        else:  # self.explicitRadio.isChecked():
            self.explicit.setEnabled(True)
            self.step.setEnabled(False)
            self.number.setEnabled(False)
            self.explicitLabel.setEnabled(True)

    def weight_show_hide(self):
        if not self.weightCheckBox.isChecked():
            self.weightLabel.show()
            self.weight.show()
        else:  # self.weightCheckBox.isChecked()
            self.weightLabel.hide()
            self.weight.hide()

    def tls_show_hide(self):
        tls_widgets = [self.tlsinFileName, self.tlsinButton, self.tls_ncyc,
                       self.tls_ncycLabel, self.tlsin_keep]
        if self.tlsCheckBox.isChecked():
            for tls_widget in tls_widgets:
                tls_widget.show()
        else:
            for tls_widget in tls_widgets:
                tls_widget.hide()

    def complete_show_hide(self):
        if self.completeCheckBox.isChecked():
            self.prerefinementCheckBox.setChecked(True)
            self.prerefinement_show_hide  # Show
            self.prerefinementCheckBox.setEnabled(False)
        else: 
            # not hide
            self.prerefinementCheckBox.setEnabled(True)
                
    def prerefinement_show_hide(self):
        prerefinement_widgets = [
            self.prerefinement_ncyc, self.prerefinement_ncycLabel, 
            self.prerefinement_addCheckBox, self.prerefinement_add,
            self.prerefinement_setCheckBox, self.prerefinement_set,
            self.prerefinement_resetCheckBox,
            self.prerefinement_shakeCheckBox, self.prerefinement_shake]
        if self.prerefinementCheckBox.isChecked():
            for prerefinement_widget in prerefinement_widgets:
                prerefinement_widget.show()
            self.prerefinement_resetCheckBox.setChecked(True)
            self.prerefinement_shakeCheckBox.setChecked(True)
        else:  # not self.prerefinementCheckBox.isChecked():
            for prerefinement_widget in prerefinement_widgets:
                prerefinement_widget.hide()

    def prerefinement_setUncheck(self):
        if self.prerefinement_resetCheckBox.isChecked():
            self.prerefinement_setCheckBox.setChecked(False)

    def prerefinement_resetUncheck(self):
        if self.prerefinement_setCheckBox.isChecked():
            self.prerefinement_resetCheckBox.setChecked(False)

    def run(self):
        # Check validity of required arguments and working directory
        if not os.path.isfile(str(self.xyzinFileName.text())) or \
                not os.path.isfile(str(self.hklinFileName.text())) or \
                not os.path.isdir(str(self.workdirName.text())):
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Structure model, merged diffracion data and "
                        "working directory must be valid!")
            # msg.setInformativeText("This is additional information")
            msg.setWindowTitle("PAIREF Error")
            msg.exec_()
            return
        # Prepare arguments
        argms = ["--XYZIN", str(self.xyzinFileName.text())]
        argms += ["--HKLIN", str(self.hklinFileName.text())]
        if self.hklin_unmergedFileName.text():
            argms += ["-u", str(self.hklin_unmergedFileName.text())]
        if self.libinFileName.text():
            argms += ["--LIBIN", str(self.libinFileName.text())]
        if self.cominFileName.text() and self.refmacRadio.isChecked():
            argms += ["--comfile", str(self.cominFileName.text())]
        elif self.defFileName.text() and self.phenixRadio.isChecked():
            argms += ["--def", str(self.defFileName.text())]
        if self.tlsinFileName.text() and self.refmacRadio.isChecked():
            argms += ["--TLSIN", str(self.tlsinFileName.text())]
        if self.refmacRadio.isChecked():
            argms += ["--refmac"]
        elif self.phenixRadio.isChecked():
            argms += ["--phenix"]
        if self.res_init.text() :
            argms += ["-i", str(self.res_init.text())]
        if self.step_numberRadio.isChecked():
            if self.step.text() :
                argms += ["-s", str(self.step.text())]
            if self.number.text():
                argms += ["-n", str(self.number.text())]
        else:  # self.explicitRadio is Checked:
            argms += ["-r", str(self.explicit.text())]
        if self.projectName.text():
            argms += ["-p", str(self.projectName.text())]
        if self.ncyc.value():
            argms += ["--ncyc", str(self.ncyc.value())]
        if self.free.text():
            argms += ["--flag", str(self.free.text())]
        if not self.weightCheckBox.isChecked() and self.refmacRadio.isChecked():
            # e.i. not auto weight
            argms += ["--weight", str(self.weight.text())]
        if self.tls_ncyc.value() and self.tlsCheckBox.isChecked() and \
                self.refmacRadio.isChecked():
            argms += ["--TLS-ncyc", str(self.tls_ncyc.value())]
        if self.tlsin_keep.isChecked() and self.tlsCheckBox.isChecked() and \
                self.refmacRadio.isChecked():
            argms += ["--TLSIN-keep"]
        if self.completeCheckBox.isChecked():
            argms += ["--complete"]
        if self.prerefinementCheckBox.isChecked():
            argms += ["--prerefinement-ncyc", str(self.prerefinement_ncyc.text())]
        if self.prerefinement_shakeCheckBox.isChecked() and self.prerefinementCheckBox.isChecked():
            argms += ["--prerefinement-shake-sites", str(self.prerefinement_shake.text())]
        if self.prerefinement_resetCheckBox.isChecked() and self.prerefinementCheckBox.isChecked():
            argms += ["--prerefinement-reset-bfactor"]
        if self.prerefinement_setCheckBox.isChecked() and self.prerefinementCheckBox.isChecked():
            argms += ["--prerefinement-set-bfactor", str(self.prerefinement_set.text())]
        if self.prerefinement_addCheckBox.isChecked() and self.prerefinementCheckBox.isChecked():
            argms += ["--prerefinement-add-to-bfactor", str(self.prerefinement_add.text())]
        if self.prerefinementCheckBox.isChecked() and \
                not self.prerefinement_shakeCheckBox.isChecked() and \
                not self.prerefinement_resetCheckBox.isChecked() and \
                not self.prerefinement_setCheckBox.isChecked() and \
                not self.prerefinement_addCheckBox.isChecked():
            argms += ["--prerefinement-no-modification"]
        # Predict new working directory
        project = str(self.projectName.text())
        if project == "":
            project = "project"
        workdir = "pairef_" + project
        while os.path.isdir(workdir):
            workdir += "_new"
        self.logfilepath = os.path.abspath(workdir) + \
            os.path.sep + "PAIREF_" + project + ".html"
        # pairef.run_pairef(argms)
        # import subprocess
        # subprocess.call([sys.executable, "-m", "pairef"] + argms)
        # Run
        argms = ["-m", "pairef"] + argms
        process = QProcess()
        print(argms)
        self.re, self.pid = process.startDetached(sys.executable, argms, ".")
        # except SystemExit:
        #     print("ignoring SystemExit")
        # Show/hide/enable buttons
        self.browserLabel.setText(
            "Job was executed, check current status in log file\n" + self.logfilepath)
        self.stopButton.show()
        self.stopButton.setEnabled(True)
        self.browserLabel.show()
        self.browserButton.show()

    def openBrowser(self):
        webbrowser.open("file://" + self.logfilepath)

    def stop(self):
        os.kill(self.pid, signal.SIGTERM)
        self.sender().setEnabled(False)
        message = self.browserLabel.text()
        message.replace("executed", "terminated")
        self.browserLabel.setText(message)


def gui():
    """Opens the window of PAIREF graphical interface (in PyQt4)"""
    app = QApplication(sys.argv)
    w = MyWindow()
    app.exec_()
