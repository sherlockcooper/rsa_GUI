#include "rsa_GUI.h"

rsa_GUI::rsa_GUI(QWidget *parent)
    : QMainWindow(parent)
{
    ui.setupUi(this);

    rsaTool = new rsa(768);
    connect(ui.generateButton, SIGNAL(clicked()), this, SLOT(generateKey()));
    connect(ui.encryptButton, SIGNAL(clicked()), this, SLOT(encryptText()));
    connect(ui.decryptButton, SIGNAL(clicked()), this, SLOT(decryptCypher()));
}

void rsa_GUI::generateKey()
{
    // ui.keyLenEdit->setText("hi");
    string content = ui.keyLenEdit->text().toStdString();
    if (content.empty())
        return;
    int len = stoi(content);
    rsaTool->setKeyLen(len);
    rsaTool->keygen();
    ostringstream helper;
    helper << "p: " << rsaTool->p << endl;
    helper << "q: " << rsaTool->q << endl;
    helper << "n: " << rsaTool->n << endl;
    helper << "phi: " << rsaTool->phi << endl;
    helper << "public key e is: " << rsaTool->e << endl;
    helper << "private key d is: " << rsaTool->d << endl;
    // helper << "time consumed: " << rsaTool->timeConsumed.count() << " ms" << endl;

    string to_show = helper.str();
    ui.keyinfoEdit->setText(QString::fromStdString(to_show));
}

void rsa_GUI::encryptText()
{
    string text = ui.textEdit->toPlainText().toStdString();
    string to_show;
    if (!rsaTool->initialized)
    {
        to_show = "rsa key not generated yet";
        
    }
    else
    {
        rsaTool->encryptStr(text, to_show);
    }
    ui.cypherEdit->setPlainText(QString::fromStdString(to_show));
    
}

void rsa_GUI::decryptCypher()
{
    string cypher = ui.cypherEdit->toPlainText().toStdString();
    string to_show;
    if (!rsaTool->initialized)
    {
        to_show = "rsa key not generated yet";

    }
    else
    {
        rsaTool->decryptStr(cypher, to_show);
    }
    ui.textEdit->setPlainText(QString::fromStdString(to_show));
}
