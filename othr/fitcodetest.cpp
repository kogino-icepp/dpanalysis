//スプリアス検出コード並びにピークサーチのコードが実データ内でしっかり機能しているかどうかを、テストシグナルを用いて検証する
#include <iostream>
#include "../headers/fitter.h"
#include <random>
using namespace std;

void fitcodetest(){
    /*
    要件定義(本関数内でやりたいこと)
    1. 実データからランダムかつ十分な量のサンプルを取り出す(ここの有意性はしっかり検証)
    2. 取り出した区間に対してピークのサイズを変えながらノイズやシグナルを載せ、それぞれ棘検出とピーク検出がどのくらいできるか実験
    */
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());
    uniform_int_distribution<> randband(1,24);//取り出すバンドをランダムに決定
    uniform_int_distribution<> randj(0,7);//取り出す測定番号を決定
    uniform_int_distribution<> randbin(sb,fb);//どの区間を引っ張り出すかを決定
}