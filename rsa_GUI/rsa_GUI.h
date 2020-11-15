#pragma once

#include <QtWidgets/QMainWindow>
#include <qpushbutton.h>
#include <qlineedit.h>
#include <qplaintextedit.h>
#include <qtextedit.h>
#include "ui_rsa_GUI.h"
#include "rsa.h"

class rsa_GUI : public QMainWindow
{
    Q_OBJECT

public:
    rsa_GUI(QWidget *parent = Q_NULLPTR);
    rsa* rsaTool;
    QPushButton* generateButton;
    QPushButton* encryptButton;
    QPushButton* decryptButton;
    QLineEdit* keyLenEdit;


private:
    Ui::rsa_GUIClass ui;

private slots:
    void generateKey();
    void encryptText();
    void decryptCypher();
};
