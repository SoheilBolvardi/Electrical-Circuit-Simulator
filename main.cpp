#include <bits/stdc++.h>
#include <fstream>

using namespace std;

class Node;

class Element;

class Node {
private:
    string name;
    double voltage;
    vector<Element *> connected_elements;
public:
    Node(string name_) : name(name_) {}

    vector<Element *> getConnectedElements() { return connected_elements; }

    string getName() { return name; }

    void addConnectedElement(Element *e) {
        connected_elements.push_back(e);
    }

    void removeConnectedElement(Element *e) {
        auto x = find(connected_elements.begin(), connected_elements.end(), e);
        if (x != connected_elements.end())
            connected_elements.erase(x);
    }

    void setName(string new_name) {
        name = new_name;
    }
};

class Circuit {
protected:
    vector<Element *> elements;
    vector<Node *> nodes;
    map<string, Node *> node_access;
    Node *gnd=nullptr;
public:
    vector<Element *> &getElements() { return elements; }

    map<string, Node *> getNodeAccess() { return node_access; }

    Node *getCreateNode(string name) {
        if (!node_access.count(name)) {
            Node *newnode = new Node(name);
            nodes.push_back(newnode);
            node_access[name] = newnode;
        }
        return node_access[name];
    }

    void addElement(Element *element) {
        elements.push_back(element);
    }

    Node* getGround() {return gnd;}

    bool setGround(string name){
        if(!node_access.count(name))
            return false;
        gnd=node_access[name];
        return true;
    }

    bool removeGround(string name){
        if(gnd!=node_access[name])
            return false;
        gnd=nullptr;
        return true;
    }

    bool renameNode(string old_name, string new_name, string& err) {
        if (node_access.find(old_name) == node_access.end()) {
            err = "ERROR: Node " + old_name + " does not exist in the circuit\n";
            return false;
        }
        if (node_access.find(new_name) != node_access.end()) {
            err = "ERROR: Node name " + new_name + " already exists\n";
            return false;
        }
        Node* node = node_access[old_name];
        node->setName(new_name);
        node_access.erase(old_name);
        node_access[new_name] = node;
        err = "SUCCESS: Node renamed from " + old_name + " to " + new_name + "\n";
        return true;
    }

};

class Element {
protected:
    string name;
    double value;
    Node *node1;
    Node *node2;
public:
    Element(string name_, double value_, Node *n1, Node *n2) :
            name(name_), value(value_), node1(n1), node2(n2) {}

    string getName() {
        return name;
    }

    virtual string getType() = 0;

    double getValue() { return value; }

    Node *getFirstNode() { return node1; }

    Node *getSecondNode() { return node2; }

    virtual ~Element() = default;
};

class Resistor : public Element {
public:
    Resistor(string name_, double value_, Node *n1, Node *n2)
            : Element(name_, value_, n1, n2) {}

    string getType() override { return "Resistor"; }
};

class Capacitor : public Element {
public:
    Capacitor(string name_, double value_, Node *n1, Node *n2)
            : Element(name_, value_, n1, n2) {}

    string getType() override { return "Capacitor"; }
};

class Inductor : public Element {
public:
    Inductor(string name_, double value_, Node *n1, Node *n2)
            : Element(name_, value_, n1, n2) {}

    string getType() override { return "Inductor"; }
};

class Diode : public Element {
protected:
    string model;
public:
    Diode(string name, Node *n1, Node *n2) :
            Element(name, 0, n1, n2), model("D") {}

    string getType() override { return "Diode"; }

    virtual double calculateCurrent(double voltage) {

    }
};

class ZenerDiode : public Diode {
private:
    double breakdownvoltage;
public:
    ZenerDiode(string name, Node *n1, Node *n2, double breakdownvoltage_ = 5)
            : Diode(name, n1, n2), breakdownvoltage(breakdownvoltage_) {
        model = "Z";
    }

    string getType() override { return "Zener Diode"; }

    double calculateCurrent(double voltage) override {

    }
};

class VoltageSource : public Element {
public:
    VoltageSource(string name_, double value_, Node *n1, Node *n2)
            : Element(name_, value_, n1, n2) {}

    virtual double getVoltage(double time) = 0;

    string getType() override { return "Voltage Source"; }
};

class DCVoltageSource : public VoltageSource {
public:
    DCVoltageSource(string name_, double value_, Node *n1, Node *n2)
            : VoltageSource(name_, value_, n1, n2) {}

    double getVoltage(double time) override {
        return getValue();
    }

    string getType() override { return "DC Voltage Source"; }
};

class SinusoidalVoltageSource : public VoltageSource {
private:
    double frequency;
    double offset;

public:
    SinusoidalVoltageSource(string name_, double amplitude, double frequency_, double offset_, Node *n1, Node *n2)
            : VoltageSource(name_, amplitude, n1, n2), frequency(frequency_), offset(offset_) {}

    double getVoltage(double time) override {
        return getValue() * sin(2 * M_PI * frequency * time) + offset;
    }

    string getType() override { return "Sinusoidal Voltage Source"; }
};

class PulseVoltageSource : public VoltageSource {
private:
    double V1;
    double V2;
    double TD;
    double TR;
    double TF;
    double TOn;
    double period;

public:
    PulseVoltageSource(string name_, double V1_, double V2_, double TD_, double TR_, double TF_, double TOn_, double period_, Node *n1, Node *n2)
            : VoltageSource(name_, V1_, n1, n2), V1(V1_), V2(V2_), TD(TD_), TR(TR_), TF(TF_), TOn(TOn_), period(period_) {}

    double getVoltage(double time) override {
        double t = fmod(time, period);

        if (t < TD) {
            return V1;
        }
        else if (t < TD + TR) {
            return V1 + (V2 - V1) * (t - TD) / TR;
        }
        else if (t < TD + TR + TOn) {
            return V2;
        }
        else if (t < TD + TR + TOn + TF) {
            return V2 - (V2 - V1) * (t - TD - TOn - TR) / TF;
        }
        else {
            return V1;
        }
    }

    string getType() override { return "Pulse Voltage Source"; }
};

class CurrentSource : public Element {
public:
    CurrentSource(string name_, double value_, Node *n1, Node *n2)
            : Element(name_, value_, n1, n2) {}

    virtual double getCurrent(double time) = 0;

    string getType() override { return "Current Source"; }
};

class DCCurrentSource : public CurrentSource {
public:
    DCCurrentSource(string name_, double value_, Node *n1, Node *n2)
            : CurrentSource(name_, value_, n1, n2) {}

    double getCurrent(double time) override {
        return getValue();
    }

    string getType() override { return "DC Current Source"; }
};

class SinusoidalCurrentSource : public CurrentSource {
private:
    double frequency;
    double offset;

public:
    SinusoidalCurrentSource(string name_, double amplitude, double frequency_, double offset_, Node *n1, Node *n2)
            : CurrentSource(name_, amplitude, n1, n2), frequency(frequency_), offset(offset_) {}

    double getCurrent(double time) override {
        return getValue() * sin(2 * M_PI * frequency * time) + offset;
    }

    string getType() override { return "Sinusoidal Current Source"; }
};

class PulseCurrentSource : public CurrentSource {
private:
    double I1;
    double I2;
    double TD;
    double TR;
    double TF;
    double TOn;
    double period;

public:
    PulseCurrentSource(string name_, double I1_, double I2_, double TD_, double TR_, double TF_, double TOn_, double period_, Node *n1, Node *n2)
            : CurrentSource(name_, I1_, n1, n2), I1(I1_), I2(I2_), TD(TD_), TR(TR_), TF(TF_), TOn(TOn_), period(period_) {}

    double getCurrent(double time) override {
        double t = fmod(time, period);

        if (t < TD) {
            return I1;
        }
        else if (t < TD + TR) {
            return I1 + (I2 - I1) * (t - TD) / TR;
        }
        else if (t < TD + TR + TOn) {
            return I2;
        }
        else if (t < TD + TR + TOn + TF) {
            return I2 - (I2 - I1) * (t - TD - TOn - TR) / TF;
        }
        else {
            return I1;
        }
    }

    string getType() override { return "Pulse Current Source"; }
};


class Controller {
public:
    string handleError(string name, Circuit *circuit, double value) {
        string type;
        switch (name[0]) {
            case 'R':
                type = "Resistor";
                break;
            case 'C':
                type = "Capacitor";
                break;
            case 'L':
                type = "Inductor";
                break;
            case 'D':
                type = "Diode";
                break;
        }
        if (value <= 0 && (type == "Resistor" || type == "Capacitor" || type == "Inductor"))
            return "Error: " + type + " cannot be zero or negative\n";
        for (auto element: circuit->getElements()) {
            if (!element)
                continue;
            if (element->getName() == name) {
                switch (name[0]) {
                    case 'R':
                        return "Resistor " + name + " already exists in the circuit\n";
                    case 'C':
                        return "Capacitor " + name + " already exists in the circuit\n";
                    case 'L':
                        return "Inductor " + name + " already exists in the circuit\n";
                    case 'D':
                        return "Diode" + name + " already exists in the circuit\n";
                }
            }
        }
        return "";
    }

    string addNewElement(string node1, string node2, string name, double value, Circuit *circuit) {
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        switch (name[0]) {
            case 'R':
                element = new Resistor(name, value, n1, n2);
                break;
            case 'L':
                element = new Inductor(name, value, n1, n2);
                break;
            case 'C':
                element = new Capacitor(name, value, n1, n2);
                break;
        }
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string removeElement(string name, Circuit *circuit) {
        auto &elements = circuit->getElements();
        for (auto it = elements.begin(); it != elements.end(); it++) {
            if ((*it)->getName() == name) {
                (*it)->getFirstNode()->removeConnectedElement(*it);
                (*it)->getSecondNode()->removeConnectedElement(*it);
                Element *deleting = *it;
                elements.erase(it);
                delete deleting;
                return name + " removed successfully!\n";
            }
        }
        switch (name[0]) {
            case 'R':
                return "Error: Cannot delete resistor; component not found\n";
            case 'L':
                return "Error: Cannot delete inductor; component not found\n";
            case 'C':
                return "Error: Cannot delete capacitor; component not found\n";
            case 'D':
                return "Error: Cannot delete diode; component not found\n";
        }
    }

    string addDiode(string node1, string node2, string name, string model, Circuit *circuit) {
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        if (model == "Z") {
            ZenerDiode *zenerDiode = new ZenerDiode(name, n1, n2);
            circuit->addElement(zenerDiode);
            n1->addConnectedElement(zenerDiode);
            n2->addConnectedElement(zenerDiode);
            return "Zener Diode added successfully!\n";
        } else if (model == "D") {
            Diode *diode = new Diode(name, n1, n2);
            circuit->addElement(diode);
            n1->addConnectedElement(diode);
            n2->addConnectedElement(diode);
            return "Diode added successfully!\n";
        }
    }

    string addGround(string node, Circuit* circuit){
        if(circuit->getGround())
            return "Error! ground node already exists in the circuit\n";
        Node *n1 = circuit->getCreateNode(node);
        circuit->setGround(node);
        return "Ground added to node " + node +" successfully!\n";
    }

    string deleteGround(string node, Circuit* circuit){
        if(!circuit->getNodeAccess().count(node))
            return "Node does not exist\n";
        if(!circuit->removeGround(node))
            return "Error! Node is not set as ground\n";
        return "Ground removed from node "+node+" successfully!\n";
    }

    void showNodes(Circuit* circuit){
        cout<<"Available nodes:\n";
        if(!circuit->getNodeAccess().size())
            cout<<"No node exists in the circuit\n";
        int index=1;
        for(auto node: circuit->getNodeAccess()){
            cout<<node.first;
            if(index!=circuit->getNodeAccess().size())
                cout<<", ";
            else
                cout<<"\n";
            index++;
        }
    }

    void list(Circuit* circuit){
        cout<<"Elements:\n";
        if(!circuit->getElements().size())
            cout<<"No elements exists in the circuit\n";
        for(auto element: circuit->getElements()){
            cout<<element->getType()<<": "<<element->getName()<<", value: "<<element->getValue()<<endl;
        }
    }

    void listElement(string type, Circuit* circuit){
        switch(type[0]){
            case 'R':
                type = "Resistor";
                break;
            case 'C':
                type = "Capacitor";
                break;
            case 'L':
                type = "Inductor";
                break;
            case 'D':
                type = "Diode";
                break;
        }
        cout<<type<<"(s):\n";
        int index=0;
        for(auto element: circuit->getElements()){
            if(element->getType()==type) {
                cout << element->getType() << ": " << element->getName() << ", value: " << element->getValue() << endl;
                index++;
            }
        }if(!index)
            cout<<"No "<<type<<" exist in the circuit\n";
    }

    string addDCVoltageSource(string name, string node1, string node2, double value, Circuit *circuit){
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        element = new DCVoltageSource(name, value, n1, n2);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addSinusoidalVoltageSource(string name, string node1, string node2, double amplitude,
                                      double frequency, double offset, Circuit *circuit){
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        element = new SinusoidalVoltageSource(name, amplitude, frequency, offset, n1, n2);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addPulseVoltageSource(string name, string node1, string node2, double V1, double V2, double TD, double TR,
                                 double TF, double TOn, double period, Circuit* circuit){
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        element = new PulseVoltageSource(name, V1, V2, TD, TR, TF, TOn, period, n1, n2);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addDCCurrentSource(string name, string node1, string node2, double value, Circuit *circuit){
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        element = new DCCurrentSource(name, value, n1, n2);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addSinusoidalCurrentSource(string name, string node1, string node2, double amplitude,
                                      double frequency, double offset, Circuit *circuit){
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        element = new SinusoidalCurrentSource(name, amplitude, frequency, offset, n1, n2);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addPulseCurrentSource(string name, string node1, string node2, double I1, double I2, double TD, double TR,
                                 double TF, double TOn, double period, Circuit* circuit){
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        element = new PulseCurrentSource(name, I1, I2, TD, TR, TF, TOn, period, n1, n2);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    void showCircuitDetails(Circuit *circuit) {
        cout << "Circuit Details:\nElements:\n";
        for (auto element: circuit->getElements()) {
            cout << "type: " << element->getType() << "   name: " << element->getName() << "  value: "
                 << element->getValue()
                 << "  node1: " << element->getFirstNode()->getName() << "  node2: "
                 << element->getSecondNode()->getName() << endl;
        }
        cout << "Node Details:\n";
        for (auto node: circuit->getNodeAccess()) {
            cout << node.first << " : " << endl;
            if (node.second->getConnectedElements().size() == 0)
                cout << "No element is connected to this node!\n";
            else {
                for (auto element: node.second->getConnectedElements())
                    cout << element->getName() << endl;
            }
        }
        cout<<"GND:\n";
        if(circuit->getGround())
            cout<<circuit->getGround()->getName()<<endl;
        else
            cout<<"No ground is specified for this circuit\n";
    }




};

class View {
private:
    Circuit *circuit;
    Controller controller;
public:
    void run() {
        circuit = new Circuit();
        string input;
        smatch match;
        regex add_element(R"(^\s*add\s+([A-Za-z])(\w+)\s+(\S+)\s+(\S+)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s*$)");
        regex remove_element(R"(^\s*delete\s+([A-Za-z])(\w+)\s*$)");
        regex add_diode(R"(^\s*add\s+(D)(\w+)\s+(\S+)\s+(\S+)\s+(D|Z)\s*$)");
        regex show_details(R"(^\s*show\s+details\s*$)");
        regex add_ground(R"(^\s*add\s+(\w+)\s+(\w+)\s*$)");
        regex delete_ground(R"(^\s*delete\s+(\w+)\s+(\w+)\s*$)");
        regex show_nodes(R"(^\s*.nodes\s*$)");
        regex list(R"(^\s*.list\s*$)");
        regex list_element(R"(^\s*.list\s+(\w+)\s*$)");
        regex rename_node(R"(^\s*\.rename\s+node\s+(\S+)\s+(\S+)\s*$)");
        regex new_file(R"(NewFile ([a-zA-Z0-9\\-_\\.:\\/()\\s]+))");
        regex addVoltageSource(R"(add\s+VoltageSource(\S+)\s+(\S+)\s+(\S+)\s+(-?\d*\.?\d+)\s*)");
        regex exit(R"(^exit$)");
        while (true) {
            getline(cin, input);
            if (regex_match(input, match, add_element)) {
                if (match[1] != "R" && match[1] != "L" && match[1] != "C") {
                    cout << "Element " << match[1] << " not found in library\n";
                    continue;
                }
                string name = match[1].str() + match[2].str();
                string node1 = match[3].str();
                string node2 = match[4].str();
                string number = match[5].str();
                string unit = match[6].str();
                double value;
                try {
                    value = stod(number);
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value\n";
                    continue;
                }
                string error = controller.handleError(name, circuit, value);
                if (error != "") {
                    cout << error;
                    continue;
                }
                if (!unit.empty()) {
                    switch (unit[0]) {
                        case 'G':
                            value *= 1e9;
                            break;
                        case 'M':
                            value *= 1e6;
                            break;
                        case 'k':
                        case 'K':
                            value *= 1e3;
                            break;
                        case 'm':
                            value *= 1e-3;
                            break;
                        case 'u':
                            value *= 1e-6;
                            break;
                        case 'n':
                            value *= 1e-9;
                            break;
                    }
                }
                cout << controller.addNewElement(node1, node2, name, value, circuit);

            } else if (regex_match(input, match, remove_element)) {
                if (match[1] != "R" && match[1] != "L" && match[1] != "C" && match[1]!="D") {
                    cout << "Element " << match[1] << " not found in library\n";
                    continue;
                }
                cout << controller.removeElement(match[1].str() + match[2].str(), circuit);
            } else if (regex_match(input, match, add_diode)) {
                string name = match[1].str() + match[2].str();
                string node1= match[3];
                string node2=match[4];
                string model=match[5];
                string error = controller.handleError(name, circuit, 0);
                if (error != "") {
                    cout << error;
                    continue;
                }
                if(match[1]!="D") {
                    cout << "Error: Element " << name << " not found in library\n";
                    continue;
                }if(model=="D" || model=="Z")
                    cout<<controller.addDiode(node1, node2, name, model, circuit);
            } else if (regex_match(input, match, add_ground)) {
                if(match[1]!="GND"){
                    cout<<"Error: Element "<<match[1]<<" not found in library\n";
                    continue;
                }
                cout << controller.addGround(match[2], circuit);
            } else if (regex_match(input, match, delete_ground)) {
                if(match[1]!="GND"){
                    cout<<"Error: Element "<<match[1]<<" not found in library\n";
                    continue;
                }
                cout << controller.deleteGround(match[2], circuit);
            } else if (regex_match(input, match, show_nodes)) {
                controller.showNodes(circuit);
            } else if (regex_match(input, match, list)) {
                controller.list(circuit);
            } else if (regex_match(input, match, list_element)) {
                if(match[1]!="R" && match[1]!="C" && match[1]!="L" && match[1]!="D" && match[1]!="Z") {
                    cout << "Error: Element " << match[1] << " not found in library\n";
                    continue;
                }controller.listElement(match[1], circuit);
            } else if (regex_match(input, match, rename_node)) {
                string old_name = match[1];
                string new_name = match[2];
                string err;
                circuit->renameNode(old_name, new_name, err);
                cout << err;
            } else if (regex_match(input, match, show_details)) {
                controller.showCircuitDetails(circuit);
            }
            else if (regex_match(input, match, new_file)) {
                string address = match[1];
                ifstream fin(address, ios::in);

                if (!fin) {
                    cerr << "Error opening file!" << endl;
                    return;
                }

                cout << "lets start" << endl;

                string line;
                regex word_regex("\\S+");
                while (getline(fin, line)) {
                    auto words_begin = sregex_iterator(line.begin(), line.end(), word_regex);
                    auto words_end = sregex_iterator();

                    vector<string> words;

                    for (auto it = words_begin; it != words_end; ++it) {
                        words.push_back(it->str());
                    }


                    string type = words[0];
                    string name = words[1];
                    string node1 = words[2];
                    string node2 = words[3];
                    string tmpValue = words[4];

                    if (node1 == "GND"){
                        controller.addGround(node1, circuit);
                    }
                    if (node2 == "GND"){
                        controller.addGround(node2, circuit);
                    }

                    double value = 0;

                    if (tmpValue[0] == '-') {
                        cout << "Error: Invalid numeric value for " << name << endl;
                        continue;
                    }

                    try {
                        if (tmpValue.back() == 'G') {
                            value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e9;
                        } else if (tmpValue.back() == 'M') {
                            value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e6;
                        }
                        else if (tmpValue.back() == 'k' || tmpValue.back() == 'K') {
                            value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e3;
                        } else if (tmpValue.back() == 'u') {
                            value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-6;
                        } else if (tmpValue.back() == 'n') {
                            value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-9;
                        } else if (tmpValue.back() == 'm') {
                            value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-3;
                        } else {
                            value = stod(tmpValue);
                        }
                    } catch (const invalid_argument &e) {
                        cout << "Error: Invalid numeric value for " << name << endl;
                        continue;
                    }

                    if (type == "R" || type == "C" || type == "L") {
                        cout << controller.addNewElement(node1, node2, name, value, circuit);

                    } else if (type == "D" || type == "Z") {
                        string model;
                        if (type == "D"){
                            model = "D";
                        }
                        if (type == "Z"){
                            model = "Z";
                        }
                        cout << controller.addDiode(node1, node2, name, model, circuit);

                    } else {
                        cout << "Element " << type << " not found in library\n";
                    }
                }

                cout << "reading file ended :)" << endl;
                fin.close();
            }


            else if (regex_match(input, match, exit)) {
                cout << "Bye Bye!\n";
                return;
            } else
                cout << "Syntax error\n";
        }
    }
};

int main() {
    View view;
    view.run();

    return 0;
}

