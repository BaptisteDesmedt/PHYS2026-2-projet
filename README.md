# 🖥️ Guide Git/GitHub pour Windows

## 📥 Installation de Git sur Windows

1. **Télécharger Git** :  
   [https://git-scm.com/download/win](https://git-scm.com/download/win)

2. **Exécuter l'installateur** :  
   - Cochez toutes les options par défaut  
   - Choisissez **Git Bash** comme terminal principal  

3. **Vérifier l'installation** :  
   Ouvrez `Git Bash` et tapez :  
   ```bash
   git --version
   # Doit afficher : git version 2.x.x.windows.1

🔧 Configuration Initiale
bash
Copy

# Définir votre identité
git config --global user.name "Votre Nom"
git config --global user.email "votre.email@domaine.com"

# Générer une clé SSH (facultatif mais recommandé)
ssh-keygen -t ed25519 -C "votre.email@domaine.com"
# Copiez le contenu de ~/.ssh/id_ed25519.pub et ajoutez-le à :
# GitHub > Settings > SSH and GPG keys

🚀 Commandes de Base pour Windows
Initialiser un projet
bash
Copy

# Créer un dossier et l'initialiser
mkdir mon-projet
cd mon-projet
git init

Cloner un dépôt existant
bash
Copy

git clone https://github.com/utilisateur/depot.git

Gérer les modifications
bash
Copy

# Ajouter tous les fichiers modifiés
git add .

# Faire un commit
git commit -m "Ajout fonctionnalité X"

# Envoyer vers GitHub
git push origin main

🛠 Workflow Typique (Git Bash)

    Créer une branche :
    bash
    Copy

    git checkout -b ma-nouvelle-branche

    Travailler et valider :
    bash
    Copy

    # Après modifications...
    git add fichier.py
    git commit -m "Correction bug login"

    Pousser vers GitHub :
    bash
    Copy

    git push origin ma-nouvelle-branche

    Ouvrir une Pull Request sur GitHub via l'interface web.

🔄 Synchroniser avec un dépôt existant
bash
Copy

# Ajouter le dépôt "officiel" comme upstream
git remote add upstream https://github.com/projet-officiel/depot.git

# Récupérer les modifications
git fetch upstream

# Fusionner avec votre branche
git merge upstream/main

⚠️ Dépannage Courant
Erreur d'authentification
bash
Copy

# Générer un token GitHub : Settings > Developer settings > Personal access tokens
git remote set-url origin https://<VOTRE_TOKEN>@github.com/utilisateur/depot.git

Annuler un fichier modifié
bash
Copy

git restore fichier.txt

Réinitialiser un dépôt
bash
Copy

git reset --hard HEAD

💡 Bonnes Pratiques Windows

    Éviter les problèmes de chemins :
    Utilisez des chemins courts sans espaces (C:/dev/mon-projet)

    Gestion des fins de ligne :
    bash
    Copy

    git config --global core.autocrlf true

    Utiliser .gitignore :
    Créez un fichier .gitignore pour exclure :
    Copy

    # Exemple pour Python
    __pycache__/
    *.pyc
    .env

📚 Ressources Utiles

    Git pour Windows Documentation

    GitHub Desktop (alternative graphique)

    Formation Git Microsoft
