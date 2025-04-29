#pragma once

#include <libAnimaCoreExport.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <array>
#include <filesystem> // C++17 ou supérieur
#include <stdexcept>
#include <limits>
#include <cmath>
#include <cstdint>
#include <cstring> // Pour std::memcpy
#include <algorithm> // Pour std::equal

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

namespace anima
{
/**
 * @class TckReader
 * @brief Permet de lire les fichiers de suivi de fibres (.tck) générés par MRtrix.
 *
 * Cette classe analyse l'en-tête ASCII et lit les données binaires des streamlines.
 * Elle supporte actuellement les types de données Float32LE et Float64LE.
 * Les streamlines sont stockées comme des vecteurs de points 3D.
 *
 * @note Utilise les fonctionnalités C++17/C++20 (notamment std::filesystem).
 * Assurez-vous de compiler avec une norme C++ appropriée (ex: -std=c++20).
 */
class TckReader
{
public:
    /**
     * @typedef Point
     * @brief Représente un point 3D (x, y, z). Utilise float par défaut.
     * Le type réel (float ou double) dépend du champ 'datatype' dans l'en-tête TCK.
     * Pour simplifier l'interface, nous utilisons double ici pour accommoder les deux,
     * mais la lecture interne respecte le type du fichier.
     * Une alternative serait de faire une classe template ou d'utiliser std::variant.
     */
    using Point = std::array<double, 3>;

    /**
     * @typedef Streamline
     * @brief Représente une streamline comme une séquence de points 3D.
     */
    using Streamline = std::vector<Point>;

    /**
     * @brief Constructeur de la classe TckReader.
     *
     * Ouvre le fichier TCK spécifié, analyse l'en-tête et lit toutes les streamlines.
     * Lève une exception std::runtime_error en cas d'échec (fichier non trouvé, format invalide, etc.).
     *
     * @param filename Chemin vers le fichier .tck.
     * @throws std::runtime_error Si le fichier ne peut pas être ouvert, lu, ou si le format est incorrect.
     */
    explicit TckReader(const std::filesystem::path& filename);

    /**
     * @brief Destructeur par défaut.
     */
    ~TckReader() = default;

    // Supprimer le constructeur par copie et l'opérateur d'affectation par copie
    TckReader(const TckReader&) = delete;
    TckReader& operator=(const TckReader&) = delete;

    // Autoriser le constructeur de déplacement et l'opérateur d'affectation par déplacement
    TckReader(TckReader&&) = default;
    TckReader& operator=(TckReader&&) = default;

    /**
     * @brief Récupère l'en-tête du fichier TCK.
     * @return Une map constante contenant les paires clé-valeur de l'en-tête.
     */
    const std::map<std::string, std::string>& getHeader() const noexcept;
    
    /**
     * @brief Récupère le nombre total de streamlines lues.
     * @return Le nombre de streamlines.
     */
	vtkPolyData *GetOutputData() {return m_OutputData;}

    /**
     * @brief Met à jour les données de sortie.
     *
     * Cette méthode est appelée pour mettre à jour les données de sortie
     * après la lecture des streamlines. Elle peut être utilisée pour
     * effectuer des opérations supplémentaires sur les données lues.
     */
    void Update();

    /**
 * @brief Convertit les données des streamlines chargées en un objet vtkPolyData.
 *
 * Chaque streamline est représentée comme une vtkPolyLine dans le vtkPolyData résultant.
 * Les points de toutes les streamlines sont stockés dans l'objet vtkPoints du vtkPolyData.
 * L'objet retourné est géré par un vtkSmartPointer pour une gestion automatique de la mémoire.
 *
 * @note Cette méthode nécessite que la bibliothèque VTK soit correctement configurée et liée
 * à votre projet lors de la compilation.
 * @return Un vtkSmartPointer pointant vers le nouvel objet vtkPolyData contenant les streamlines.
 * Retourne un PolyData vide (mais valide) si aucune streamline n'a été chargée.
 * @warning Par défaut, vtkPoints utilise des 'float' pour stocker les coordonnées.
 * Si vos données TCK sont en 'double' (Float64LE), une conversion (et potentielle perte
 * de précision) aura lieu lors de l'insertion des points. Pour conserver la précision
 * double, vous pouvez décommenter la ligne `points->SetDataType(VTK_DOUBLE);`
 * dans l'implémentation de cette méthode.
 */
    vtkSmartPointer<vtkPolyData> convertToVtkPolyData() const;

protected:
    /**
     * @brief Récupère toutes les streamlines lues depuis le fichier.
     * @return Un vecteur constant contenant toutes les streamlines. Chaque streamline est un vecteur de points.
     * @warning Pour les fichiers TCK très volumineux, cette méthode peut consommer beaucoup de mémoire.
     */
    const std::vector<Streamline>& getStreamlines() const noexcept;

    /**
     * @brief Récupère le chemin du fichier TCK lu.
     * @return Le chemin du fichier.
     */
    const std::filesystem::path& getFilePath() const noexcept;

    /**
     * @brief Récupère la valeur d'une clé spécifique de l'en-tête.
     * @param key La clé à rechercher dans l'en-tête.
     * @return Une chaîne contenant la valeur associée à la clé.
     * @throws std::out_of_range si la clé n'est pas trouvée dans l'en-tête.
     */
    const std::string& getHeaderValue(const std::string& key) const;


private:
    /**
     * @enum DataType
     * @brief Énumère les types de données supportés pour les coordonnées.
     */
    enum class SupportedDataType
	{
        FLOAT32_LE, /**< Simple précision (4 octets), Little Endian */
        FLOAT64_LE, /**< Double précision (8 octets), Little Endian */
        UNSUPPORTED /**< Type de données non supporté */
    };

    /**
     * @brief Analyse l'en-tête ASCII du fichier TCK.
     *
     * Lit les lignes jusqu'à trouver le marqueur "END". Valide le format,
     * extrait les paires clé-valeur et détermine le type de données et le décalage
     * vers les données binaires.
     *
     * @param file Le flux de fichier ouvert en mode binaire.
     * @throws std::runtime_error Si l'en-tête est mal formé ou si le type de données n'est pas supporté.
     */
    void parseHeader(std::ifstream& file);

    /**
     * @brief Lit les données binaires des streamlines depuis le fichier.
     *
     * Se positionne au début des données binaires (déterminé par parseHeader)
     * et lit les coordonnées jusqu'à la fin du fichier, en les regroupant
     * par streamline.
     *
     * @param file Le flux de fichier ouvert en mode binaire.
     * @throws std::runtime_error En cas d'erreur de lecture ou de fin de fichier inattendue.
     * @throws std::bad_alloc Si l'allocation mémoire échoue (pour les très gros fichiers).
     */
    void readStreamlinesData(std::ifstream& file);

    /**
     * @brief Lit une streamline unique en fonction du type de données.
     * @tparam T Le type de données des coordonnées (float ou double).
     * @param file Le flux de fichier binaire.
     * @param currentStreamline Le vecteur où stocker les points de la streamline lue.
     * @return true si une streamline complète a été lue (terminée par NaN/Inf), false si EOF est atteint avant la fin.
     */
    template <typename T>
    bool readSingleStreamline(std::ifstream& file, Streamline& currentStreamline);

    /**
     * @brief Vérifie si un triplet de coordonnées représente la fin d'une streamline.
     *
     * La fin est marquée par (NaN, NaN, NaN) ou (Inf, Inf, Inf).
     * @param coords Le triplet de coordonnées à vérifier.
     * @return true si le triplet marque la fin, false sinon.
     */
    static bool isEndOfStreamline(const std::array<double, 3>& coords) noexcept;
    static bool isEndOfStreamline(const std::array<float, 3>& coords) noexcept;


    std::filesystem::path m_filename;                 ///< Chemin du fichier TCK.
    std::map<std::string, std::string> m_header;      ///< Contenu de l'en-tête (clé-valeur).
    std::vector<Streamline> m_streamlines;            ///< Stockage de toutes les streamlines.
    std::streampos m_dataOffset = 0;                  ///< Position dans le fichier où les données binaires commencent.
    SupportedDataType m_dataType = SupportedDataType::UNSUPPORTED; ///< Type de données détecté.
    size_t m_pointSize = 0;                           ///< Taille en octets d'un point (3 * sizeof(type)).
    size_t m_coordSize = 0;                           ///< Taille en octets d'une coordonnée (sizeof(type)).
	

private:
    vtkSmartPointer <vtkPolyData> m_OutputData;
    std::string m_FileName;
};
}