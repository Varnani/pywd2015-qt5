from PyQt5 import QtWidgets, QtGui
from gui import configurespots_widget
from src import constants
from configparser import ConfigParser
from src.helpers import methods, messenger


class Widget(QtWidgets.QWidget, configurespots_widget.Ui_SpotConfigureWidget):
    def __init__(self, parent):
        super(Widget, self).__init__()
        self.setupUi(self)

        self.setWindowIcon(QtGui.QIcon(constants.MAIN_ICON_PATH))

        self.main_window = parent

        self.star1_treewidget = SpotTreeWidget(self.star1_groupbox, self, 1)
        self.star2_treewidget = SpotTreeWidget(self.star2_groupbox, self, 2)

        self.a_btn_group = QtWidgets.QButtonGroup()
        self.b_btn_group = QtWidgets.QButtonGroup()

        self.connect_signals()

    def connect_signals(self):
        self.whatsthis_btn.clicked.connect(QtWidgets.QWhatsThis.enterWhatsThisMode)
        self.clear_a_b_button.clicked.connect(self.clear_a_b)
        self.spotconfigsave_btn.clicked.connect(self.save_spots)
        self.spotconfigload_btn.clicked.connect(self.load_spots)

    def clear_a_b(self):
        self.a_btn_group.setExclusive(False)
        self.b_btn_group.setExclusive(False)

        for btn in self.a_btn_group.buttons():
            btn.setChecked(False)

        for btn in self.b_btn_group.buttons():
            btn.setChecked(False)

        self.a_btn_group.setExclusive(True)
        self.b_btn_group.setExclusive(True)

    def write_into_parser(self, parser):
        parser.add_section(constants.CONFIG_SECTION_SPOT_PARAMS)
        parser.set(constants.CONFIG_SECTION_SPOT_PARAMS, "fspot1", str(self.fspot1_ipt.value()))
        parser.set(constants.CONFIG_SECTION_SPOT_PARAMS, "fspot2", str(self.fspot2_ipt.value()))
        parser.set(constants.CONFIG_SECTION_SPOT_PARAMS, "ifsvm1", str(self.ifsmv1_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_SPOT_PARAMS, "ifsvm2", str(self.ifsmv2_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_SPOT_PARAMS, "kspev", str(self.kspev_groupbox.isChecked()))
        parser.set(constants.CONFIG_SECTION_SPOT_PARAMS, "nomax", str(self.nomax_combobox.currentIndex()))
        parser.set(constants.CONFIG_SECTION_SPOT_PARAMS, "kspot", str(self.kspot_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_SPOT_PARAMS, constants.CONFIG_SECTION_STAR1_SPOTS + "count",
                   str(self.star1_treewidget.tree_widget.invisibleRootItem().childCount()))
        parser.set(constants.CONFIG_SECTION_SPOT_PARAMS, constants.CONFIG_SECTION_STAR2_SPOTS + "count",
                   str(self.star2_treewidget.tree_widget.invisibleRootItem().childCount()))

        parser = self.star1_treewidget.write_into_parser(parser)
        parser = self.star2_treewidget.write_into_parser(parser)

        return parser

    def read_from_parser(self, parser):
        self.fspot1_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SPOT_PARAMS, "fspot1"))
        self.fspot2_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SPOT_PARAMS, "fspot2"))
        self.ifsmv1_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_SPOT_PARAMS, "ifsvm1"))
        self.ifsmv2_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_SPOT_PARAMS, "ifsvm2"))
        self.kspev_groupbox.setChecked(parser.getboolean(constants.CONFIG_SECTION_SPOT_PARAMS, "kspev",))
        self.nomax_combobox.setCurrentIndex(parser.getint(constants.CONFIG_SECTION_SPOT_PARAMS, "nomax",))
        self.kspot_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_SPOT_PARAMS, "kspot",))

        self.star1_treewidget.read_from_parser(parser)
        self.star2_treewidget.read_from_parser(parser)

    def save_spots(self):
        save_path = methods.save_file(self, suffix=".pywdspot", name_filter="PyWD2015 spot file (*.pywdspot)")
        if save_path is not None:
            parser = self.write_into_parser(ConfigParser())

            with open(save_path, "w") as destination:
                parser.write(destination)

            msg = messenger.Messenger("info", "Save completed:")
            msg.set_info(save_path)
            msg.show()

    def load_spots(self):
        load_path = methods.load_file(self, suffix=".pywdspot", name_filter="PyWD2015 spot file (*.pywdspot);;"
                                                                            "PyWD2015 project file (*.pywdproject)")
        if load_path is not None:
            parser = ConfigParser()
            with open(load_path, "r") as source:
                parser.read_file(source)

            self.read_from_parser(parser)

            msg = messenger.Messenger("info", "Load completed:")
            msg.set_info(load_path)
            msg.show()

    def get_all_spots(self, get_ab=False):
        def _from_item_to_list(item, parent, ab=False):
            value_list = [parent.itemWidget(item, 2).value(), parent.itemWidget(item, 3).value(),
                          parent.itemWidget(item, 4).value(), parent.itemWidget(item, 5).value(),
                          parent.itemWidget(item, 6).value(), parent.itemWidget(item, 7).value(),
                          parent.itemWidget(item, 8).value(), parent.itemWidget(item, 9).value()]

            if ab:
                value_list.append(parent.itemWidget(item, 0).isChecked())
                value_list.append(parent.itemWidget(item, 1).isChecked())

            return value_list

        s1_spot_items = self.star1_treewidget.get_all_items()
        s2_spot_items = self.star2_treewidget.get_all_items()

        s1_spots = [_from_item_to_list(spot_item, self.star1_treewidget.tree_widget, ab=get_ab)
                    for spot_item in s1_spot_items]
        s2_spots = [_from_item_to_list(spot_item, self.star2_treewidget.tree_widget, ab=get_ab)
                    for spot_item in s2_spot_items]

        return s1_spots, s2_spots


class SpotTreeWidget:
    def __init__(self, groupbox, parent, star):

        self.star = star
        self.groupbox = groupbox
        self.v_layout = QtWidgets.QVBoxLayout(self.groupbox)
        self.tree_widget = QtWidgets.QTreeWidget()
        self.v_layout.addWidget(self.tree_widget)
        self.h_layout = QtWidgets.QHBoxLayout()
        self.add_spot_button = QtWidgets.QPushButton()
        self.add_spot_button.setText("Add Spot")
        self.remove_spot_button = QtWidgets.QPushButton()
        self.remove_spot_button.setText("Remove Selected Spot")
        self.h_layout.addWidget(self.add_spot_button)
        self.h_layout.addWidget(self.remove_spot_button)
        self.v_layout.addLayout(self.h_layout)

        self.parent_widget = parent

        columns = ["A", "B", "LAT", "LON", "RADSP", "TEMSP", "Tstart", "Tmax1", "Tmax2", "Tend"]

        header_item = QtWidgets.QTreeWidgetItem()
        for index, column in enumerate(columns):
            header_item.setText(index, column)

        self.tree_widget.setHeaderItem(header_item)
        self.header_item_count = len(columns)

        self.tree_widget.header().setSectionResizeMode(1)
        self.tree_widget.header().setSectionsMovable(False)

        self.connect_signals()

    def connect_signals(self):
        self.add_spot_button.clicked.connect(self.insert_row)
        self.remove_spot_button.clicked.connect(self.remove_spot_row)

    def remove_spot_row(self):
        item = self.selected_item()
        self.tree_widget.headerItem().removeChild(item)

    def selected_item(self):
        selecteditem = self.tree_widget.selectedItems()
        if len(selecteditem) > 0:
            return selecteditem[0]
        else:
            return None

    def insert_row(self, ):
        item = QtWidgets.QTreeWidgetItem(self.tree_widget)

        a_radio_button = QtWidgets.QRadioButton()
        a_radio_button.setText("A")
        b_radio_button = QtWidgets.QRadioButton()
        b_radio_button.setText("B")

        self.parent_widget.a_btn_group.addButton(a_radio_button)
        self.parent_widget.b_btn_group.addButton(b_radio_button)

        self.tree_widget.setItemWidget(item, 0, a_radio_button)
        self.tree_widget.setItemWidget(item, 1, b_radio_button)

        i = 2
        while i < self.header_item_count:
            spinbox = QtWidgets.QDoubleSpinBox(self.tree_widget)
            spinbox.setButtonSymbols(2)
            spinbox.setDecimals(8)
            spinbox.setMaximum(9999999)
            spinbox.setMinimum(0)
            self.tree_widget.setItemWidget(item, i, spinbox)
            i = i + 1

        return item

    def get_section(self):

        section = ""

        if self.star == 1:
            section = constants.CONFIG_SECTION_STAR1_SPOTS

        if self.star == 2:
            section = constants.CONFIG_SECTION_STAR2_SPOTS

        return section

    def get_all_items(self):
        count = self.tree_widget.invisibleRootItem().childCount()
        children = []
        i = 0
        while i < count:
            children.append(self.tree_widget.invisibleRootItem().child(i))
            i = i + 1

        return children

    def write_into_parser(self, parser):
        section = self.get_section()

        for i, child in enumerate(self.get_all_items()):
            sct = section + str(i + 1)
            parser.add_section(sct)
            parser.set(sct, "a", str(self.tree_widget.itemWidget(child, 0).isChecked()))
            parser.set(sct, "b", str(self.tree_widget.itemWidget(child, 1).isChecked()))
            parser.set(sct, "lat", str(self.tree_widget.itemWidget(child, 2).value()))
            parser.set(sct, "lon", str(self.tree_widget.itemWidget(child, 3).value()))
            parser.set(sct, "radsp", str(self.tree_widget.itemWidget(child, 4).value()))
            parser.set(sct, "temsp", str(self.tree_widget.itemWidget(child, 5).value()))
            parser.set(sct, "tstart", str(self.tree_widget.itemWidget(child, 6).value()))
            parser.set(sct, "tmax1", str(self.tree_widget.itemWidget(child, 7).value()))
            parser.set(sct, "tmax2", str(self.tree_widget.itemWidget(child, 8).value()))
            parser.set(sct, "tend", str(self.tree_widget.itemWidget(child, 9).value()))

        return parser

    def read_from_parser(self, parser):
        self.tree_widget.clear()

        i = 1
        while i <= parser.getint(constants.CONFIG_SECTION_SPOT_PARAMS, self.get_section() + "count"):
            section = self.get_section() + str(i)
            item = self.insert_row()
            self.tree_widget.itemWidget(item, 0).setChecked(parser.getboolean(section, "a"))
            self.tree_widget.itemWidget(item, 1).setChecked(parser.getboolean(section, "b"))
            self.tree_widget.itemWidget(item, 2).setValue(parser.getfloat(section, "lat"))
            self.tree_widget.itemWidget(item, 3).setValue(parser.getfloat(section, "lon"))
            self.tree_widget.itemWidget(item, 4).setValue(parser.getfloat(section, "radsp"))
            self.tree_widget.itemWidget(item, 5).setValue(parser.getfloat(section, "temsp"))
            self.tree_widget.itemWidget(item, 6).setValue(parser.getfloat(section, "tstart"))
            self.tree_widget.itemWidget(item, 7).setValue(parser.getfloat(section, "tmax1"))
            self.tree_widget.itemWidget(item, 8).setValue(parser.getfloat(section, "tmax2"))
            self.tree_widget.itemWidget(item, 9).setValue(parser.getfloat(section, "tend"))

            i = i + 1
